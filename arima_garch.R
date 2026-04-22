
library(rugarch)

setwd("C:/Users/julie/Desktop/Trondheim S7/Time Series/project/")

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
as_num <- function(x) as.numeric(x)

# ---------- Load & prepare ----------
data <- read.csv("boe_fx.csv")
data <- na.omit(data[, c("XUDLNKD", "XUDLUSS")])
Serie <- as_num(data$XUDLNKD / data$XUDLUSS)

n_test <- 200
train_lev <- head(Serie, -n_test)  # levels
test_lev  <- tail(Serie,  n_test)  # levels

stopifnot(all(is.finite(train_lev)), all(is.finite(test_lev)))

# ---------- Helpers ----------
# back-transform a one-step forecast from d-fold differences to level
invert_one_step <- function(last_levels, d, pred_diff1) {
  # last_levels: tail levels needed (vector of length >= d+1 ideally)
  # pred_diff1: forecast of Δ^d x_{t+1}
  L <- length(last_levels)
  x_t   <- last_levels[L]
  if (d == 0) {
    # model is on levels; pred is already the level forecast
    return(pred_diff1)
  } else if (d == 1) {
    # x_{t+1} = x_t + Δ x_{t+1}
    return(x_t + pred_diff1)
  } else if (d == 2) {
    # x_{t+1} = 2 x_t - x_{t-1} + Δ^2 x_{t+1}
    x_tm1 <- last_levels[L - 1]
    return(2 * x_t - x_tm1 + pred_diff1)
  } else if (d == 3) {
    # x_{t+1} = 3 x_t - 3 x_{t-1} + x_{t-2} + Δ^3 x_{t+1}
    x_tm1 <- last_levels[L - 1]
    x_tm2 <- last_levels[L - 2]
    return(3 * x_t - 3 * x_tm1 + x_tm2 + pred_diff1)
  } else {
    stop("d > 3 not supported in this utility")
  }
}

# Build fixed.pars with exactly the names the spec expects
build_fixed_pars <- function(coefs, pm, qm, pv, qv, include_mean) {
  expected <- c(
    if (include_mean) "mu",
    # variance params
    "omega",
    if (pv > 0) paste0("alpha", seq_len(pv)),
    if (qv > 0) paste0("beta",  seq_len(qv)),
    # mean ARMA params
    if (pm > 0) paste0("ar", seq_len(pm)),
    if (qm > 0) paste0("ma", seq_len(qm))
  )
  keep <- intersect(names(coefs), expected)
  as.list(coefs[keep])
}

cat_line <- function(...) cat(sprintf(...), "\n")

# ---------- Grid search ----------
# Mean side: d in 0..3, ARMA(p_mean, q_mean)
# Variance side: sGARCH(p_var, q_var)
d_grid     <- 0:3
pm_grid    <- 0:7
qm_grid    <- 0:7
pv_grid    <- 0:3
qv_grid    <- 0:3

# --- Fit ARIMA and rolling forecast ---
min_rmse_err    <- Inf
min_aic_err    <- Inf
best_rmse_pdq_pq  <- c(NA_integer_, NA_integer_, NA_integer_,NA_integer_, NA_integer_)
best_aic_pdq_pq  <- c(NA_integer_, NA_integer_, NA_integer_,NA_integer_, NA_integer_)
cat_err    <- function(...) cat(sprintf(...), "\n")

results <- data.frame(
  da = integer(),
  pa = integer(),
  qa = integer(),
  pg = integer(),
  qg = integer(),
  RMSE = numeric(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)
for (da in d_grid) {
  # differenced series used for *fitting*
  # we’ll always forecast levels by back-transforming
  train_y <- if (da == 0) train_lev else diff(train_lev, differences = da)
  if (length(train_y) < 30) next  # too short, skip
  stopifnot(all(is.finite(train_y)))
  
  
  for (pa in pm_grid) {
    for (qa in qm_grid) {
      for (pg in pv_grid) {
        for (qg in qv_grid) {

          if (pg == 0 && qg == 0) next
          if (pa == 0 && qa == 0) next
          
          # 1) Fit once (estimate params) on *train_y*
          base_spec <- ugarchspec(
            variance.model = list(model = "sGARCH", garchOrder = c(pg, qg)),
            mean.model     = list(armaOrder = c(pa, qa), include.mean = TRUE, arfima = FALSE),
            distribution.model = "norm"
          )
          
          fit <- tryCatch(
            ugarchfit(spec = base_spec, data = train_y, solver.control = list(trace = 0)),
            error = function(e) { cat_line("  -> fit failed: %s", e$message); NULL }
          )
          if (is.null(fit)) next
          
          # 2) Freeze parameters
          coefs <- coef(fit)
          inc_mu <- "mu" %in% names(coefs)
          fixed_pars <- build_fixed_pars(coefs, pa, qa, pg, qg, inc_mu)
          
          fixed_spec <- ugarchspec(
            variance.model = list(model = "sGARCH", garchOrder = c(pg, qg)),
            mean.model     = list(armaOrder = c(pa, qa), include.mean = inc_mu, arfima = FALSE),
            distribution.model = fit@model$modeldesc$distribution,
            fixed.pars     = fixed_pars
          )
          
          # 3) Rolling one-step forecasts on *levels*, via differencing & back-transform
          preds <- numeric(length(test_lev))
          ok <- TRUE
          
          # We reveal test_lev progressively on the *level* scale
          for (i in seq_along(test_lev)) {
            # history in levels up to time t = end of train + (i-1) test points
            lev_hist <- if (i == 1) train_lev else c(train_lev, test_lev[seq_len(i - 1)])
            
            # create differenced series of order d for the model input
            y_hist <- if (da == 0) lev_hist else diff(lev_hist, differences = da)
            
            if (!length(y_hist) || any(!is.finite(y_hist))) {
              cat_line("  -> bad hist at i=%d", i); ok <- FALSE; break
            }
            
            fc <- tryCatch(
              ugarchforecast(fitORspec = fixed_spec, data = y_hist, n.ahead = 1),
              error = function(e) { cat_line("  -> forecast failed at i=%d: %s", i, e$message); NULL }
            )
            if (is.null(fc)) { ok <- FALSE; break }
            
            # forecast on the differenced scale (or level if d=0)
            pred_diff1 <- as_num(fc@forecast$seriesFor[1])
            
            # back-transform to level
            # need last d+1 levels to invert
            need <- max(1, da + 1L)
            if (length(lev_hist) < need) { ok <- FALSE; cat_line("  -> not enough levels to invert at i=%d", i); break }
            
            tail_lev <- tail(lev_hist, need)
            preds[i] <- invert_one_step(tail_lev, da, pred_diff1)
          }
          
          if (!ok) { cat_line("  -> skip: rolling step broke"); next }
          
          # 4) Score
          rmse_err <- rmse(test_lev, preds)
          aic_value = infocriteria(fit)[1]
          
          cat_err("ARIMA(p=%d,d=%d,q=%d)-GARCH(p=%d,q=%d) RMSE = %.6f AIC = %.6f",
                  pa, da, qa,pg,qg, rmse_err,aic_value) 
          
          if (rmse_err < min_rmse_err) {
            min_rmse_err   <- rmse_err
            best_rmse_pdq_pq <- c(pa, da, qa,pg,qg)
            cat_err("Better rmse")
          }
          if (aic_value < min_aic_err){
            min_aic_err <- aic_value
            best_aic_pdq_pq <- c(pa,da,qa,pg,qg)
            cat_err("Better aic")
          }
          
          results <- rbind(
            results,
            data.frame(
              da = da,
              pa = pa,
              qa = qa,
              pg = pg,
              qg = qg,
              RMSE = rmse_err,
              AIC = aic_value
            )
          )
          
        }
      }
    }
  }
}

cat_err("The best ARIMA-GARCH model for RMSE error is ARIMA(p=%d,d=%d,q=%d)-GARCH(p=%d,q=%d) with RMSE = %.6f",
        best_rmse_pdq_pq[1],best_rmse_pdq_pq[2],best_rmse_pdq_pq[3],
        best_rmse_pdq_pq[4],best_rmse_pdq_pq[5],  min_rmse_err)
cat_err("The best ARIMA-GARCH model for AIC error is ARIMA(p=%d,d=%d,q=%d)-GARCH(p=%d,q=%d) with AIC = %.6f",
        best_aic_pdq_pq[1],best_aic_pdq_pq[2],best_aic_pdq_pq[3],
        best_aic_pdq_pq[4],best_aic_pdq_pq[5],  min_aic_err)

write.csv(results, "ARIMA-GARCH-results.csv", row.names = TRUE)