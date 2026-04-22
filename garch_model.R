library(rugarch)

setwd("C:/Users/julie/Desktop/Trondheim S7/Time Series/project/")

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
cat_line <- function(...) cat(sprintf(...), "\n")

# --- Load and prepare data ---
data <- read.csv("boe_fx.csv")
data <- na.omit(data[, c("XUDLNKD", "XUDLUSS")])

# Level series to predict (ratio)
Serie <- data$XUDLNKD / data$XUDLUSS

# Train/test split on LEVELS
n_test   <- 200
train_lev <- head(Serie, -n_test)
test_lev  <- tail(Serie,  n_test)

# --- PRECOMPUTE RETURNS ONCE (log-returns) ---
logret <- function(x) diff(log(as.numeric(x)))
r_all     <- logret(Serie)                                 # length = length(Serie)-1
n_train_r <- length(train_lev) - 1                         # returns count in train
train_ret <- r_all[1:n_train_r]
test_ret  <- r_all[(n_train_r + 1):length(r_all)]          # aligned to test_lev (same length as test_lev)

stopifnot(length(test_ret) == length(test_lev))

# --- Grid search over GARCH variance orders (p,q) on returns ---
min_rmse_err    <- Inf
min_aic_err    <- Inf
best_rmse_pq  <- c(NA_integer_, NA_integer_)
best_aic_pq  <- c(NA_integer_, NA_integer_)
cat_err    <- function(...) cat(sprintf(...), "\n")

results <- data.frame(
  p = integer(),
  q = integer(),
  RMSE = numeric(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

for (p in 0:5) {
  for (q in 0:5) {
    if (p==0 && q==0) next
    
    # 1) Spec and fit once on TRAIN RETURNS
    spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
      mean.model     = list(armaOrder = c(0,0), include.mean = TRUE),
      distribution.model = "norm"
    )
    
    fit <- tryCatch(
      ugarchfit(spec = spec, data = train_ret, solver.control = list(trace = 0)),
      error = function(e) NULL
    )
    if (is.null(fit)) { cat_line("  -> skip: fit failed"); next }
    
    # 2) Freeze parameters (no re-estimation later)
    pars <- coef(fit)
    spec_fixed <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
      mean.model     = list(armaOrder = c(0,0), include.mean = TRUE),
      distribution.model = fit@model$modeldesc$distribution,
      fixed.pars     = as.list(pars)
    )
    
    preds_lev <- numeric(length(test_lev))
    ok <- TRUE
    
    # 3) Rolling one-step MEAN forecast on RETURNS, then invert to LEVEL
    # At step i, history in LEVELS is: train_lev + first (i-1) test_lev points
    for (i in seq_along(test_lev)) {
      
      # build returns history by slicing, no recompute from levels
      y_hist_ret <- if (i == 1) train_ret else c(train_ret, test_ret[seq_len(i - 1)])
      if (!length(y_hist_ret) || any(!is.finite(y_hist_ret))) {
        ok <- FALSE; cat_line("  -> bad returns history at i=%d", i); break
      }
      
      # Deterministic 1-step forecast (conditional mean) on returns, params frozen
      fc <- tryCatch(
        ugarchforecast(fitORspec = spec_fixed, data = y_hist_ret, n.ahead = 1),
        error = function(e) NULL
      )
      if (is.null(fc)) { ok <- FALSE; cat_line("  -> forecast failed at i=%d", i); break }
      
      # Predicted next return (mean of the predictive distribution)
      r_hat <- as.numeric(fc@forecast$seriesFor[1])
      
      # Invert to level: x_{t+1} = x_t * exp(r_{t+1})
      x_last <- if (i == 1) tail(train_lev, 1) else test_lev[i - 1]
      x_hat  <- x_last * exp(r_hat)
      
      preds_lev[i] <- x_hat
    }
    
    if (!ok) { cat_line("  -> skip: rolling step failed"); next }
    
    # 4) Score
    rmse_err <- rmse(test_lev, preds_lev)
    aic_value = infocriteria(fit)$AIC
    cat_err("GARCH(p=%d,q=%d) RMSE = %.6f AIC = %.6f",
            p, q, rmse_err,aic_value) 
    
    if (rmse_err < min_rmse_err) {
      min_rmse_err   <- rmse_err
      best_rmse_pq <- c(p, q)
      cat_err("Better rmse")
    }
    if (aic_value < min_aic_err){
      min_aic_err <- aic_value
      best_aic_pq <- c(p,q)
      cat_err("Better aic")
    }
    
    results <- rbind(
      results,
      data.frame(
        p = p,
        q = q,
        RMSE = rmse_err,
        AIC = aic_value
      )
    )
    
  }
}

cat_err("The best GARCH model for RMSE error is GARCH(p=%d,q=%d) with RMSE = %.6f",
        best_rmse_pq[1], best_rmse_pq[2], min_rmse_err)
cat_err("The best GARCH model for RMSE error is GARCH(p=%d,q=%d) with AIC = %.6f",
        best_aic_pq[1], best_aic_pq[2], min_aic_err)

write.csv(results, "GARCH-results.csv", row.names = TRUE)
