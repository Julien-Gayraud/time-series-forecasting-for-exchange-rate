library(forecast)
setwd("C:/Users/julie/Desktop/Trondheim S7/Time Series/project/")
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))

# --- Load and prepare data ---
data <- read.csv("boe_fx.csv")

# Drop missing rows for the two columns
data <- na.omit(data[, c("XUDLNKD", "XUDLUSS")])

# Create the ratio series
Serie <- data$XUDLNKD / data$XUDLUSS

# --- Train/test split ---
n_test <- 200
train <- head(Serie, -n_test)
test  <- tail(Serie,  n_test)

# --- Fit ARIMA and rolling forecast ---
min_rmse_err    <- Inf
min_aic_err    <- Inf
best_rmse_pdq  <- c(NA_integer_, NA_integer_, NA_integer_)
best_aic_pdq  <- c(NA_integer_, NA_integer_, NA_integer_)
cat_err    <- function(...) cat(sprintf(...), "\n")

results <- data.frame(
  d = integer(),
  p = integer(),
  q = integer(),
  RMSE = numeric(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

for (d in 0:3) {
  for (p in 0:7) {
    for (q in 0:7) {
      # 1) Fit once to get parameters
      fit <- tryCatch(
        Arima(train, order = c(p, d, q), method = "ML"),
        error = function(e) NULL
      )
      
      if (is.null(fit)) { cat_err("  -> skip: fit failed"); next }
      
      
      # 2) Freeze parameters (no refit later), keep same mean/drift setting
      par_fixed <- coef(fit)
      inc_mean  <- isTRUE(fit$include.mean) || "intercept" %in% names(par_fixed)
      
      state_model <- tryCatch(
        Arima(train,
              order = c(p, d, q),
              fixed = par_fixed,
              include.mean = inc_mean,
              method = "ML",
              transform.pars = FALSE),
        error = function(e) NULL
      )
      if (is.null(state_model)) { cat_err("  -> skip: state init failed"); next }
      
      # 3) Rolling one-step-ahead forecast over test (update state only)
      preds <- numeric(length(test))
      ok <- TRUE
      for (i in seq_along(test)) {
        
        fc <- tryCatch(
          as.numeric(forecast(state_model, h = 1)$mean[1]),
          error = function(e) NA_real_
        )
        
        if (is.na(fc)) { ok <- FALSE; break }
        preds[i] <- fc
        
        y_hist <- c(train, test[seq_len(i)])
        state_model <- tryCatch(
          Arima(y_hist,
                order = c(p, d, q),
                fixed = par_fixed,
                include.mean = inc_mean,
                method = "ML",
                transform.pars = FALSE),
          error = function(e) NULL
        )
        if (is.null(state_model)) { ok <- FALSE; break }
      }
      if (!ok) { cat_err("  -> skip: rolling update failed"); next }
      
      # 4) Score
      rmse_err <- rmse(test, preds)
      aic_value = AIC(fit)

      cat_err("ARIMA(p=%d,d=%d,q=%d) RMSE = %.6f AIC = %.6f",
        p, d, q, rmse_err,aic_value) 
      
      if (rmse_err < min_rmse_err) {
        min_rmse_err <- rmse_err
        best_rmse_pdq <- c(p, d, q)
        cat_err("Better rmse")
      }
      if (aic_value < min_aic_err){
        min_aic_err <- aic_value
        best_aic_pdq <- c(p,d,q)
        cat_err("Better aic")
      }
      
      results <- rbind(
        results,
        data.frame(
          d = d,
          p = p,
          q = q,
          RMSE = rmse_err,
          AIC = aic_value
        )
      )

      
      
    }
  }
}
write.csv(results, "ARIMA-results.csv", row.names = TRUE)

cat_err("The best ARIMA model for RMSE error is ARIMA(%d,%d,%d) with RMSE = %.6f",
        best_rmse_pdq[1], best_rmse_pdq[2], best_rmse_pdq[3], min_rmse_err)
cat_err("The best ARIMA model for AIC value is ARIMA(%d,%d,%d) with AIC = %.6f",
        best_aic_pdq[1], best_aic_pdq[2], best_aic_pdq[3], min_aic_err)

