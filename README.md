# Model Selection for Forecasting the USD/NOK Exchange Rate

## 1. Context

This project was carried out at **NTNU** as part of the course **Time Series**.  
The aim of the course was to apply forecasting and statistical modeling methods to real temporal datasets.

---

## 2. Project Overview

The project focuses on the **USD/NOK exchange rate**.  

The objective was to compare the performance of several models and determine which one is the most suitable for next-day forecasting.

Models tested:

- **ARIMA**
- **GARCH**
- **ARIMA-GARCH**

Evaluation metrics:

- **RMSE**: prediction accuracy  
- **AIC**: model quality / complexity trade-off

The final conclusions of the project are summarized in the poster:  
**JulienGayraud_Poster_TimeSeries.pdf**

---

## 3. Project Files

- **boe_fx.csv** → raw exchange rate dataset
- **arima_model.R** → ARIMA modeling  
- **garch_model.R** → GARCH modeling  
- **arima_garch.R** → hybrid ARIMA-GARCH modeling  
- **ARIMA-results.csv** → ARIMA results  
- **GARCH-results.csv** → GARCH results  
- **ARIMA-GARCH-results.csv** → hybrid model results  
- **ResultsAnalysis.ipynb** → final analysis and visualizations  
- **JulienGayraud_Poster_TimeSeries.pdf** → final project summary poster
