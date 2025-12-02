# Beyond GARCH: Volatility Forecasting with Copulas

This project investigates whether incorporating tail dependence via copulas improves volatility forecasting for the S&P 500.  
It contains:
- The main script with all model estimation and forecasting routines  
- A PDF report with full methodology, results, and figures :contentReference[oaicite:0]{index=0}

## Summary
We estimate and compare multiple models:
- Standard GARCH and asymmetric GARCH  
- Realized GARCH (symmetric and asymmetric)  
- GARCH-X models including VIX and EPU  
- Copula-augmented Realized GARCH (Gaussian, Student-t, Clayton, Gumbel)

Using S&P 500 daily data (2000–2024), we evaluate 1-, 5-, and 21-day volatility forecasts using MSE and QLIKE.  
Key result: While Student-t copulas capture strong tail dependence in-sample, they **do not improve out-of-sample forecasts**, and the **asymmetric Realized GARCH** model remains the most accurate.

## Files
- `script.R` — all estimation, copula fitting, and forecasting code  
- `Case_Report.pdf` — full report with tables, plots, and interpretation  
