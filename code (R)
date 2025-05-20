# ==========================================================
# In‐sample ML Estimation for SGARCH | AGARCH | SRGARCH | ARGARCH
# Hand‐coded likelihood + stats4::mle / optim to replicate Prof’s two‐step MATLAB procedure
# ==========================================================
library(stats4)
library(xts)
library(dplyr)
library(readxl)
library(MASS)   # for ginv()
library(knitr) 

# 0) Read & prepare data ---------------------------------------------
df <- read_excel("data.xlsx") %>%
  setNames(c("Date","Open","Close","LogReturn","RealizedVar","VIX","EPU")) %>%
  mutate(
    Date = as.Date(as.character(Date), "%Y%m%d"),
    r    = LogReturn,
    x    = RealizedVar
  ) %>%
  dplyr::select(Date, r, x) %>%
  na.omit()

rt <- df$r
xt <- df$x
n  <- length(rt)

# 1) SGARCH(1,1) -----------------------------------------------------
nll_sgarch <- function(mu, omega, alpha, beta) {
  if (alpha + beta >= 1) return(1e10)
  sigma2 <- rep(omega/(1-alpha-beta), n)
  for (t in 2:n) {
    eps      <- rt[t-1] - mu
    sigma2[t] <- omega + alpha * eps^2 + beta * sigma2[t-1]
  }
  -sum(dnorm(rt, mean = mu, sd = sqrt(sigma2), log = TRUE))
}
m_sg <- mle(
  nll_sgarch,
  start  = list(mu = mean(rt), omega = 0.02, alpha = 0.10, beta = 0.85),
  method = "L-BFGS-B",
  lower  = c(mu = -Inf, omega = 1e-6, alpha = 1e-6, beta = 1e-6),
  upper  = c(mu =  Inf, omega = Inf,  alpha = 1,     beta = 1),
  control= list(maxit = 1e4)
)
cat("\n--- SGARCH(1,1) Estimates ---\n")
print(summary(m_sg))

# 2) AGARCH(1,1) Glosten‐GJR via optim ----------------------------
nll_agarch <- function(par) {
  mu     <- par[1]
  omega  <- par[2]
  alpha1 <- par[3]
  alpha2 <- par[4]
  beta   <- par[5]
  if ((alpha1 + alpha2)/2 + beta >= 1) return(1e10)
  sigma2 <- rep(var(rt), n)
  for (t in 2:n) {
    eps       <- rt[t-1] - mu
    sigma2[t] <- omega +
      alpha1 * (eps^2) * (eps < 0) +
      alpha2 * (eps^2) * (eps >= 0) +
      beta   * sigma2[t-1]
    if (!is.finite(sigma2[t]) || sigma2[t] <= 0) return(1e10)
  }
  -sum(dnorm(rt, mean = mu, sd = sqrt(sigma2), log = TRUE))
}
start_ag <- c(mu = mean(rt), omega = 0.02, alpha1 = 0.05, alpha2 = 0.05, beta = 0.90)
low_ag   <- c(mu = -Inf, omega = 1e-6, alpha1 = 1e-6, alpha2 = 1e-6, beta = 1e-6)
up_ag    <- c(mu =  Inf, omega = Inf, alpha1 = 1,     alpha2 = 1,     beta = 1)

res_ag <- optim(
  par     = start_ag,
  fn      = nll_agarch,
  method  = "L-BFGS-B",
  lower   = low_ag,
  upper   = up_ag,
  hessian = TRUE
)
vcov_ag <- tryCatch(solve(res_ag$hessian), error = function(e) ginv(res_ag$hessian))
se_ag    <- sqrt(diag(vcov_ag))
coef_ag  <- res_ag$par
logLik_ag<- -res_ag$value
cat("\n--- AGARCH(1,1) Estimates ---\n")
print(round(coef_ag, 4))
cat("Std. errors:\n"); print(round(se_ag, 4))
cat(sprintf("LogLik = %.4f\n", logLik_ag))

# 3) SRGARCH: Symmetric Realized GARCH -----------------------------
# 3.1) Partial likelihood for returns
nll_r_SRG <- function(mu, omega, alpha, beta, gamma) {
  if (alpha + beta >= 1) return(1e10)
  sig0   <- (omega + gamma * mean(xt)) / (1 - alpha + - beta)
  sigma2 <- rep(sig0, n)
  for (t in 2:n) {
    sigma2[t] <- omega +
      alpha * (rt[t-1] - mu)^2 +
      beta  * sigma2[t-1] +
      gamma * xt[t-1]
  }
  -sum(dnorm(rt, mean = mu, sd = sqrt(sigma2), log = TRUE))
}
m_r_SRG <- mle(
  nll_r_SRG,
  start  = list(mu = mean(rt), omega = 0.02, alpha = 0.10, beta = 0.85, gamma = 0.05),
  method = "L-BFGS-B",
  lower  = c(mu=-Inf,omega=1e-6,alpha=1e-6,beta=1e-6,gamma=1e-6),
  upper  = c(mu= Inf,omega=Inf,alpha=1,     beta=1,     gamma=1),
  control= list(maxit=1e4)
)
pars_r <- coef(m_r_SRG)
cat("\n--- Partial likelihood l(r) for SRGARCH ---\n"); print(pars_r)

# 3.2) Measurement equation regression xt ~ sigma^2 + (eps^2 − 1)
sigma2_r <- rep((pars_r["omega"] + pars_r["gamma"] * mean(xt)) /
                  (1 - pars_r["alpha"] - pars_r["beta"]), n)
for (t in 2:n) {
  sigma2_r[t] <- pars_r["omega"] +
    pars_r["alpha"] * (rt[t-1] - pars_r["mu"])^2 +
    pars_r["beta"]  * sigma2_r[t-1] +
    pars_r["gamma"] * xt[t-1]
}
eps_r <- (rt - pars_r["mu"]) / sqrt(sigma2_r)
dt_reg <- data.frame(xt = xt, s2 = sigma2_r, e2m1 = eps_r^2 - 1)
coef_reg <- coef(lm(xt ~ s2 + e2m1, data = dt_reg))
cat("\n--- Measurement equation (SRGARCH) ---\n"); print(coef_reg)

# 3.3) Joint likelihood l(r) + l(x|r)
nll_SRG_full <- function(mu, omega, alpha, beta, gamma, xi, phi, tau2, su2) {
  if (alpha + beta >= 1 || su2 <= 0) return(1e10)
  sig0   <- (omega + gamma * mean(xt)) / (1 - alpha - beta)
  sigma2 <- rep(sig0, n); u <- numeric(n)
  for (t in 2:n) {
    sigma2[t] <- omega +
      alpha * (rt[t-1] - mu)^2 +
      beta  * sigma2[t-1] +
      gamma * xt[t-1]
    eps      <- (rt[t] - mu) / sqrt(sigma2[t])
    u[t]     <- xt[t] - xi - phi * sigma2[t] - tau2 * (eps^2 - 1)
  }
  llr <- sum(dnorm(rt, mean = mu, sd = sqrt(sigma2), log = TRUE))
  llx <- sum(dnorm(u,   mean = 0,  sd = sqrt(su2),    log = TRUE))
  -(llr + llx)
}
start_full_SRG <- c(
  mu    = pars_r["mu"],
  omega = pars_r["omega"],
  alpha = pars_r["alpha"],
  beta  = pars_r["beta"],
  gamma = pars_r["gamma"],
  xi    = coef_reg["(Intercept)"],
  phi   = coef_reg["s2"],
  tau2  = coef_reg["e2m1"],
  su2   = var(residuals(lm(xt ~ s2 + e2m1, dt_reg)))
)
lower_full_SRG <- c(mu=-Inf,omega=1e-6,alpha=1e-6,beta=1e-6,
                    gamma=1e-6,xi=-Inf,  phi=1e-6,  tau2=1e-6, su2=1e-6)
upper_full_SRG <- c(mu= Inf,omega=Inf, alpha=1,    beta=1,
                    gamma=1,   xi= Inf,   phi=Inf,   tau2=Inf, su2=Inf)
m_full_SRG <- mle(
  nll_SRG_full,
  start  = start_full_SRG,
  method = "L-BFGS-B",
  lower  = lower_full_SRG,
  upper  = upper_full_SRG,
  control= list(maxit=1e4)
)
cat("\n--- Full SRGARCH Estimates ---\n"); print(summary(m_full_SRG))

# 4) ARGARCH: Asymmetric Realized GARCH via optim -----------------
# 4.1) Prepare sigma^2 & eps from SR step
pars_rA   <- pars_r
sigma2_rA <- sigma2_r
eps_rA    <- eps_r

# 4.2) Measurement equation regression xt ~ s2 + eps + (eps^2 - 1)
dtA  <- data.frame(xt = xt,
                   s2 = sigma2_rA,
                   eps = eps_rA,
                   e2m1 = eps_rA^2 - 1)
regA <- lm(xt ~ s2 + eps + e2m1, data = dtA); coefA <- coef(regA)
if (is.na(coefA["s2"])) {
  regA <- lm(xt ~ eps + e2m1, data = dtA); coefA <- coef(regA)
  phi0 <- 0.1
} else {
  phi0 <- coefA["s2"]
}
xi0   <- coefA["(Intercept)"]
tau10 <- ifelse(is.na(coefA["eps"]),   0, coefA["eps"])
tau20 <- ifelse(is.na(coefA["e2m1"]), 0.1, coefA["e2m1"])
su20  <- var(residuals(regA)); if (su20 <= 0) su20 <- var(xt) / 2

start_full_ARG <- c(
  mu     = pars_rA["mu"],
  omega  = pars_rA["omega"],
  alpha  = pars_rA["alpha"],
  beta   = pars_rA["beta"],
  gamma  = pars_rA["gamma"],
  xi     = xi0,
  phi    = phi0,
  tau1   = tau10,
  tau2   = tau20,
  su2    = su20
)
lower_full_ARG <- c(mu=-Inf,omega=1e-6,alpha=1e-6,beta=1e-6,
                    gamma=1e-6,xi=-Inf,phi=1e-6,
                    tau1=-Inf, tau2=1e-6, su2=1e-6)
upper_full_ARG <- c(mu= Inf,omega= Inf,alpha=1,    beta=1,
                    gamma=1,   xi= Inf,  phi=Inf,
                    tau1=Inf,  tau2=Inf,su2=Inf)

obj_full_ARG <- function(par) {
  mu     <- par[1]; omega <- par[2]; alpha <- par[3]
  beta   <- par[4]; gamma <- par[5]; xi    <- par[6]
  phi    <- par[7]; tau1  <- par[8]; tau2  <- par[9]
  su2    <- par[10]
  if (alpha+beta >= 1 || su2 <= 0) return(1e10)
  sig0   <- (omega + gamma * mean(xt)) / (1 - alpha - beta)
  sigma2 <- rep(sig0, n); u <- numeric(n)
  for (t in 2:n) {
    sigma2[t] <- omega +
      alpha * (rt[t-1] - mu)^2 +
      beta  * sigma2[t-1] +
      gamma * xt[t-1]
    eps_t     <- (rt[t] - mu) / sqrt(sigma2[t])
    u[t]      <- xt[t] - xi - phi * sigma2[t] - tau1 * eps_t - tau2 * (eps_t^2 - 1)
    if (!is.finite(sigma2[t]) || sigma2[t] <= 0) return(1e10)
  }
  llr <- sum(dnorm(rt, mean = mu, sd = sqrt(sigma2), log = TRUE))
  llx <- sum(dnorm(u,   mean = 0,  sd = sqrt(su2),    log = TRUE))
  -(llr + llx)
}

res_full_ARG <- optim(
  par     = start_full_ARG,
  fn      = obj_full_ARG,
  method  = "L-BFGS-B",
  lower   = lower_full_ARG,
  upper   = upper_full_ARG,
  hessian = TRUE,
  control = list(trace = 0, maxit = 1e4)
)
coef_full   <- res_full_ARG$par
vcov_full   <- tryCatch(solve(res_full_ARG$hessian), error = function(e) ginv(res_full_ARG$hessian))
se_full     <- sqrt(diag(vcov_full))
logLik_full <- -res_full_ARG$value

cat("\n--- Full ARGARCH Estimates (optim) ---\n")
print(round(coef_full, 4))
cat("Std. errors:\n"); print(round(se_full, 4))
cat(sprintf("LogLik = %.4f\n", logLik_full), "\n")

lr_sg   <- as.numeric(logLik(m_sg))       # SGARCH returns only
lf_sg   <- lr_sg                          # no measurement eq

lr_ag   <- logLik_ag                      # AGARCH returns only
lf_ag   <- lr_ag                          # no measurement eq

lr_srg  <- as.numeric(logLik(m_r_SRG))    # SRGARCH returns only
lf_srg  <- as.numeric(logLik(m_full_SRG)) # SRGARCH joint

lr_arg  <- lr_srg                         # ARGARCH returns only
lf_arg  <- logLik_full                    # ARGARCH joint

res_tab <- data.frame(
  Model        = c("SGARCH",  "AGARCH",   "SRGARCH",    "ARGARCH"),
  `l(r)`       = c(lr_sg,     lr_ag,      lr_srg,       lr_arg),
  `l(r)+l(x|r)` = c(lf_sg,     lf_ag,      lf_srg,       lf_arg)
)
res_tab[,2:3] <- round(res_tab[,2:3], 2)

kable(res_tab,
      col.names = c("Model", "l(r)", "l(r)+l(x|r)"),
      caption   = "Partial vs Joint Log‐Likelihood for Four Models")

# ===================================================================
# 6) PLOTS: returns, VIX, ΔVIX vs returns, σ path, variance measures
# ===================================================================
library(ggplot2)
library(scales)
library(lubridate)

df_plot <- df %>%
  mutate(
    sigma    = as.numeric(sqrt(sigma2_r)),  # use SRG sigma² for illustration; swap as needed
    diff_vix = c(NA, diff(x))               # x here is RealizedVar; if VIX needed, adjust
  )

theme_set(
  theme_classic(base_family = "Times") +
    theme(text = element_text(size = 20))
)
year_breaks <- seq(ymd("2000-01-01"), ymd("2024-12-31"), by = "1 year")

# Plot 1: log returns
ggplot(df_plot, aes(Date, r)) +
  geom_line() +
  scale_x_date(breaks = year_breaks, labels = date_format("%y")) +
  coord_cartesian(xlim = c(ymd("2000-01-01"), ymd("2024-12-31")),
                  ylim = c(min(rt), max(rt))) +
  labs(y = expression(r[t])) -> p1
print(p1)

# Plot 2: VIX (using xt if available)
ggplot(df_plot, aes(Date, x)) +  # x is RealizedVar; replace with VIX if added to df_plot
  geom_line() +
  scale_x_date(breaks = year_breaks, labels = date_format("%y")) +
  coord_cartesian(xlim = c(ymd("2000-01-01"), ymd("2024-12-31")),
                  ylim = c(0, max(df$x))) +
  labs(y = "RealizedVar") -> p2
print(p2)

# Plot 3: ΔVariance vs returns
ggplot(filter(df_plot, !is.na(diff_vix)),
       aes(x = diff_vix, y = r)) +
  geom_point() +
  labs(x = expression(Delta~x[t]), y = expression(r[t])) -> p3
print(p3)

# Plot 4: returns + sigma path
ggplot(df_plot, aes(Date)) +
  geom_line(aes(y = r)) +
  geom_line(aes(y = sigma), color = "red") +
  scale_x_date(breaks = year_breaks, labels = date_format("%y")) +
  coord_cartesian(xlim = c(ymd("2000-01-01"), ymd("2024-12-31")),
                  ylim = c(min(c(rt, df_plot$sigma)), max(c(rt, df_plot$sigma)))) +
  annotate("text", x = ymd("2005-01-01"), y = max(rt), label = "data & sigma", hjust = 0) -> p4
print(p4)

# Plot 5: variance measures
ggplot(df_plot, aes(Date)) +
  geom_point(aes(y = rt^2), shape = 1) +
  geom_line(aes(y = x)) +
  geom_line(aes(y = sigma^2), color = "red", linetype = "dotted") +
  scale_x_date(breaks = year_breaks, labels = date_format("%y")) +
  coord_cartesian(xlim = c(ymd("2000-01-01"), ymd("2024-12-31")),
                  ylim = c(0, max(c(df_plot$x, df_plot$sigma^2, rt^2)))) +
  labs(y = "Variance measures") -> p5
print(p5)























# 扩张窗口
# ── Rolling Out‐of‐Sample Forecasting for 4 Models ──────────────────────────────

# 0) Load packages
library(readxl)
library(zoo)
library(rugarch)
library(openxlsx)

# 1) Read & prepare data
df_raw <- read_excel("data.xlsx")
names(df_raw) <- c("Date","Open","Close","LogReturn","x","VIX","Policy")
df_raw$Date <- as.Date(as.character(df_raw$Date), "%Y%m%d")
# Subset and rename
df <- df_raw[, c("Date","LogReturn","x","VIX")]
names(df)[names(df)=="LogReturn"] <- "r"

# 2) Split into in‐sample and out‐of‐sample
n  <- nrow(df)
n0 <- floor(0.5 * n)
test_dates <- df$Date[(n0+1):n]

# 3) Forecast horizons
horizons <- c(1, 5, 21)
H <- length(horizons)

# Precompute HAR regressors for full sample
RV_d <- df$x
RV_w <- rollapplyr(df$x,  5, mean, fill = NA)
RV_m <- rollapplyr(df$x, 21, mean, fill = NA)

# Prepare storage
N_oos <- n - n0
PV_vix_mat  <- matrix(NA, N_oos, H, dimnames=list(NULL, paste0("VIX_d", horizons)))
PV_har_mat  <- matrix(NA, N_oos, H, dimnames=list(NULL, paste0("HAR_d", horizons)))
PV_ag_mat   <- matrix(NA, N_oos, H, dimnames=list(NULL, paste0("AGARCH_d", horizons)))
PV_arg_mat  <- matrix(NA, N_oos, H, dimnames=list(NULL, paste0("ARGARCH_d", horizons)))

# 4) AGARCH spec for rolling
spec_ag <- ugarchspec(
  variance.model     = list(model="gjrGARCH", garchOrder=c(1,1)),
  mean.model         = list(armaOrder=c(0,0)),
  distribution.model = "norm"
)

# 5) ARGARCH parameters (from your in‐sample two‐step MLE)
pars  <- coef_full
omega <- pars["omega.omega"]; alpha <- pars["alpha.alpha"]
beta  <- pars["beta.beta"];  gamma <- pars["gamma.gamma"]
mu    <- pars["mu.mu"]
xi    <- pars["xi.(Intercept)"]; phi <- pars["phi.s2"]
tau1  <- pars["tau1.eps"];    tau2 <- pars["tau2.e2m1"]

# Precompute state up to n0 for ARGARCH
sigma2_all <- numeric(n); eps_all <- numeric(n)
sigma2_all[1] <- (omega + gamma*mean(df$x[1:n0]))/(1 - alpha - beta)
eps_all[1]    <- df$r[1] - mu
for (t in 2:n0) {
  eps_all[t]    <- df$r[t] - mu
  sigma2_all[t] <- omega +
    alpha * eps_all[t-1]^2 +
    beta  * sigma2_all[t-1] +
    gamma * df$x[t-1]
}

# 6) Rolling loop
for (i in 1:N_oos) {
  t0 <- n0 + i - 1  # last in‐sample index
  
  # 6.1) VIX dynamic
  for (j in 1:H) {
    d <- horizons[j]
    PV_vix_mat[i,j] <- (d/250) * df$VIX[t0]
  }
  
  # 6.2) HAR‐RV dynamic
  for (j in 1:H) {
    d <- horizons[j]
    # build training sample for this horizon
    end_idx <- t0 - d
    idx     <- which(!is.na(RV_w[1:end_idx]) & !is.na(RV_m[1:end_idx]))
    idx     <- intersect(idx, 1:end_idx)
    TV <- sapply(idx, function(t) mean(df$x[(t+1):(t+d)]))
    reg_df <- data.frame(
      TV   = TV,
      RV_d = RV_d[idx],
      RV_w = RV_w[idx],
      RV_m = RV_m[idx]
    )
    fit <- lm(TV ~ RV_d + RV_w + RV_m, data = reg_df)
    # regressors at time t0
    PV_har_mat[i,j] <- fit$coef[1] +
      fit$coef[2] * RV_d[t0] +
      fit$coef[3] * RV_w[t0] +
      fit$coef[4] * RV_m[t0]
  }
  
  # 6.3) AGARCH rolling via manual fit + forecast
  fit_ag_i <- ugarchfit(spec_ag, data = df$r[1:t0], solver="solnp", fit.control=list(stationarity=1))
  fc_ag    <- ugarchforecast(fit_ag_i, n.ahead=max(horizons))
  sigma_fc <- as.numeric(fc_ag@forecast$sigmaFor)
  for (j in 1:H) {
    PV_ag_mat[i,j] <- sigma_fc[horizons[j]]^2
  }
  
  # 6.4) ARGARCH rolling
  s2_now <- sigma2_all[t0]; e_now <- eps_all[t0]; x_now <- df$x[t0]
  for (j in 1:H) {
    d <- horizons[j]
    s2 <- s2_now; e <- e_now; xprev <- x_now
    for (h in 1:d) {
      s2    <- omega + alpha*e^2 + beta*s2 + gamma*xprev
      e     <- 0
      xprev <- xi + phi*s2 + tau1*e + tau2*(e^2 - 1)
    }
    PV_arg_mat[i,j] <- xprev
  }
  # update ARGARCH state for next day
  eps_all[t0+1]    <- df$r[t0+1] - mu
  sigma2_all[t0+1] <- omega +
    alpha * eps_all[t0]^2 +
    beta  * sigma2_all[t0] +
    gamma * df$x[t0]
}

# 7) Assemble output and write to Excel
oos <- data.frame(Date = df$Date[(n0+1):n])
oos <- cbind(
  oos,
  as.data.frame(PV_vix_mat),
  as.data.frame(PV_har_mat),
  as.data.frame(PV_ag_mat),
  as.data.frame(PV_arg_mat)
)
names(oos) <- c(
  "Date",
  paste0("VIX_d", horizons),
  paste0("HAR_d", horizons),
  paste0("AGARCH_d", horizons),
  paste0("ARGARCH_d", horizons)
)
oos$Realized <- df$x[(n0+1):n]

wb <- createWorkbook()
addWorksheet(wb, "OOS")
writeData(wb, "OOS", oos)
saveWorkbook(wb, "Expanding Window.xlsx", overwrite = TRUE)






# ── Q5: 计算 4×3×2×2=48 个平均损失 (MSE & QLIKE) ─────────────────────────────

# 0) 加载所需包
library(readxl)
library(dplyr)
library(tidyr)
library(knitr)
library(openxlsx)

# 1) 读入“扩张窗口”结果
df <- read_excel("C:/Users/BXand/Downloads/Expanding Window_pilot1.xlsx") %>%
  mutate(Date = as.Date(Date))  # 如有必要把 Date 转成 Date 类型

# 2) 参数设定
horiz   <- c(1, 5, 21)
factor2 <- 1/0.7

# 3) 构造两类目标变量 TV1（r^2） 和 TV2（realized，乘 1/0.7）
for(d in horiz) {
  df <- df %>%
    mutate(
      !!paste0("TV1_d", d) := Reduce(`+`, lapply(1:d, function(i) lead(r^2, i))),
      !!paste0("TV2_d", d) := factor2 * Reduce(`+`, lapply(1:d, function(i) lead(Realized, i)))
    )
}

# 4) 损失函数
mse   <- function(f, h) (f - h)^2
qlike <- function(f, h) log(h) + f / h

# 5) Pivot to long, 对齐 horizon，再计算 Loss
long <- df %>%
  # 5.1 将 12 条预测拉长
  pivot_longer(
    cols         = matches("^(VIX|HAR|AGARCH|ARGARCH)_d"),
    names_to     = c("Model", "d"),
    names_pattern= "(.*)_d(\\d+)",
    values_to    = "Forecast"
  ) %>%
  mutate(d = as.integer(d)) %>%
  # 5.2 将 TV1/TV2 拉长
  pivot_longer(
    cols         = matches("^TV[12]_d"),
    names_to     = c("Target", "d2"),
    names_pattern= "(TV[12])_d(\\d+)",
    values_to    = "TargetValue"
  ) %>%
  mutate(d2 = as.integer(d2)) %>%
  # 5.3 只保留预测 horizon 与目标 horizon 对应的行
  filter(d == d2) %>%
  # 5.4 删除中间列 d2（注意这里显式调用 dplyr::select 以避免被其他包遮蔽）
  dplyr::select(-d2) %>%
  # 5.5 计算两种损失
  mutate(
    Loss_MSE   = mse(Forecast, TargetValue),
    Loss_QLIKE = qlike(Forecast, TargetValue)
  )

# 6) 汇总平均损失
summary <- long %>%
  group_by(Model, d, Target) %>%
  summarise(
    MSE   = mean(Loss_MSE,   na.rm = TRUE),
    QLIKE = mean(Loss_QLIKE, na.rm = TRUE),
    .groups = "drop"
  )

# 7) 打印漂亮的 Markdown 表格
for(tar in c("TV1", "TV2")) {
  cat("\n## Average losses for", tar, "\n")
  tab <- summary %>%
    filter(Target == tar) %>%
    arrange(Model, d) %>%
    pivot_longer(
      cols       = c(MSE, QLIKE),
      names_to   = "Loss",
      values_to  = "Value"
    ) %>%
    pivot_wider(
      names_from = c("d", "Loss"),
      values_from= "Value"
    )
  print(
    kable(
      tab,
      format  = "markdown",
      digits  = 4,
      caption = paste0("Avg. losses (", tar, ")")
    )
  )
}

# 8) （可选）把结果写到 Excel
wb <- createWorkbook()
addWorksheet(wb, "Losses_TV1")
addWorksheet(wb, "Losses_TV2")
writeData(wb, "Losses_TV1", summary %>% filter(Target=="TV1"))
writeData(wb, "Losses_TV2", summary %>% filter(Target=="TV2"))
saveWorkbook(wb, "losses_summary.xlsx", overwrite = TRUE)
































# ==========================================================
# Bonus Q6: 1-Step Ahead Variance Forecast Comparison
#   Models:
#     1. Standard GARCH(1,1)
#     2. GARCH-X (VIX)
#     3. GARCH-X (EPU)
#     4. Realized GARCH
# ==========================================================

# 0) Load required libraries ---------------------------------
library(readxl)    # read_excel()
library(rugarch)   # ugarchspec(), ugarchroll()
library(xts)       # xts()
library(ggplot2)   # autoplot()

# 1) Read & clean data ---------------------------------------
raw <- read_excel("data.xlsx")
names(raw) <- c("Date","spx_open","spx_close","r",
                "RealizedVar","VIX","EPU")

df <- data.frame(
  Date = as.Date(as.character(raw$Date), "%Y%m%d"),
  r    = raw$r,
  x    = raw$RealizedVar,
  # VIX is annualized vol in percent (e.g. 27.01)
  VIX  = raw$VIX / 100,
  EPU  = raw$EPU
)
df <- na.omit(df)

# 2) Convert to xts for rugarch ------------------------------
df_xts <- xts(df[, c("r","x","VIX","EPU")], order.by = df$Date)

# 3) Specify models ------------------------------------------
# 3.1 Standard GARCH(1,1)
spec_sgarch <- ugarchspec(
  variance.model = list(model="sGARCH", garchOrder=c(1,1)),
  mean.model     = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model="norm"
)

# 3.2 GARCH-X with VIX
spec_gx_vix <- ugarchspec(
  variance.model = list(
    model="sGARCH", garchOrder=c(1,1),
    external.regressors = as.matrix(df_xts$VIX)
  ),
  mean.model     = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model="norm"
)

# 3.3 GARCH-X with EPU
spec_gx_epu <- ugarchspec(
  variance.model = list(
    model="sGARCH", garchOrder=c(1,1),
    external.regressors = as.matrix(df_xts$EPU)
  ),
  mean.model     = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model="norm"
)

# 3.4 Realized GARCH using model="realGARCH" in variance.model
spec_realgarch <- ugarchspec(
  variance.model = list(model="realGARCH", garchOrder=c(1,1)),
  mean.model     = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model="norm"
)

# 4) Rolling one-step forecasts ------------------------------
n_oos        <- 250     # out-of-sample size
refit_every  <- 25      # refit frequency
refit_window <- "moving"

# helper to run ugarchroll; pass realizedVol for RealGARCH
run_roll <- function(spec, useRV=FALSE) {
  args <- list(
    spec            = spec,
    data            = df_xts$r,
    forecast.length = n_oos,
    refit.every     = refit_every,
    refit.window    = refit_window,
    keep.coef       = TRUE
  )
  if(useRV) args$realizedVol <- df_xts$x
  do.call(ugarchroll, args)
}

roll_sg    <- run_roll(spec_sgarch)
roll_vix   <- run_roll(spec_gx_vix)
roll_epu   <- run_roll(spec_gx_epu)
roll_real  <- run_roll(spec_realgarch, useRV=TRUE)

# 5) Extract σ² forecasts ------------------------------------
# convert roll objects to data.frames
df_sg   <- as.data.frame(roll_sg)
df_vix  <- as.data.frame(roll_vix)
df_epu  <- as.data.frame(roll_epu)
df_real <- as.data.frame(roll_real)

# retrieve 1-step ahead σ² (Sigma column is conditional sd)
f_sg_vec   <- df_sg$Sigma^2
f_vix_vec  <- df_vix$Sigma^2
f_epu_vec  <- df_epu$Sigma^2
f_real_vec <- df_real$Sigma^2

# true realized variance from original xts
realized_oos <- tail(df_xts$x, n_oos)

# time index for OOS
idx_oos <- tail(index(df_xts), n_oos)

# convert vectors into xts series
f_sg   <- xts(f_sg_vec,   order.by=idx_oos)
f_vix  <- xts(f_vix_vec,  order.by=idx_oos)
f_epu  <- xts(f_epu_vec,  order.by=idx_oos)
f_real <- xts(f_real_vec, order.by=idx_oos)

# 6) Compute MSE and QLIKE -----------------------------------
mse    <- function(f, r) mean((coredata(f) - coredata(r))^2)
qlike <- function(f, r) mean(log(coredata(f)) + coredata(r) / coredata(f))

results <- data.frame(
  Model = c("GARCH(1,1)",
            "GARCH-X (VIX)",
            "GARCH-X (EPU)",
            "RealGARCH"),
  MSE   = c(
    mse(f_sg,   realized_oos),
    mse(f_vix,  realized_oos),
    mse(f_epu,  realized_oos),
    mse(f_real, realized_oos)
  ),
  QLIKE = c(
    qlike(f_sg,   realized_oos),
    qlike(f_vix,  realized_oos),
    qlike(f_epu,  realized_oos),
    qlike(f_real, realized_oos)
  )
)
print(results)

# 7) Plot forecasts vs realized ------------------------------
# Merge all series and then drop any NA rows
plot_xts <- merge(realized_oos, f_sg, f_vix, f_epu, f_real)
plot_xts <- na.omit(plot_xts)
colnames(plot_xts) <- c("Realized", results$Model)

autoplot(plot_xts) +
  labs(
    title = "One-Day-Ahead Variance Forecasts",
    x     = "Date",
    y     = expression(hat(sigma)[t+1]^2),
    color = "Model"
  ) +
  theme_minimal()

















# ==========================================================
# Bonus Q6 (extended): Full‐Sample Conditional Variance Plots
#   Fit each model on the entire sample, extract σ²_t, and plot
#   alongside realized variance x_t over the full sample.
# ==========================================================

# 0) Load libraries ------------------------------------------
library(readxl)    # read_excel()
library(rugarch)   # ugarchspec(), ugarchfit(), sigma()
library(xts)       # xts()
library(ggplot2)   # autoplot()

# 1) Read & prepare data -------------------------------------
raw <- read_excel("data.xlsx")
names(raw) <- c("Date","spx_open","spx_close","r",
                "RealizedVar","VIX","EPU")

df <- data.frame(
  # convert YYYYMMDD → Date
  Date = as.Date(as.character(raw$Date), "%Y%m%d"),
  r    = raw$r,
  x    = raw$RealizedVar,
  # VIX: percent → decimal
  VIX  = raw$VIX / 100,
  EPU  = raw$EPU
)
df <- na.omit(df)

# convert to xts for rugarch
df_xts <- xts(df[, c("r","x","VIX","EPU")], order.by=df$Date)

# 2) Specify model specs on full sample ---------------------
# 2.1 Standard sGARCH(1,1)
spec_sgarch <- ugarchspec(
  variance.model     = list(model="sGARCH", garchOrder=c(1,1)),
  mean.model         = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model = "norm"
)

# 2.2 GARCH-X with VIX
spec_gx_vix <- ugarchspec(
  variance.model     = list(
    model="sGARCH", garchOrder=c(1,1),
    external.regressors = as.matrix(df_xts$VIX)
  ),
  mean.model         = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model = "norm"
)

# 2.3 GARCH-X with EPU
spec_gx_epu <- ugarchspec(
  variance.model     = list(
    model="sGARCH", garchOrder=c(1,1),
    external.regressors = as.matrix(df_xts$EPU)
  ),
  mean.model         = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model = "norm"
)

# 2.4 Realized GARCH (measurement equation)
spec_realgarch <- ugarchspec(
  variance.model     = list(model="realGARCH", garchOrder=c(1,1)),
  mean.model         = list(armaOrder=c(0,0), include.mean=TRUE),
  distribution.model = "norm"
)

# 3) Fit each model on the full sample -----------------------
# 3.1 Fit standard GARCH
fit_sg   <- ugarchfit(spec_sgarch,    data=df_xts$r)

# 3.2 Fit GARCH-X (VIX)
fit_vix  <- ugarchfit(spec_gx_vix,    data=df_xts$r)

# 3.3 Fit GARCH-X (EPU)
fit_epu  <- ugarchfit(spec_gx_epu,    data=df_xts$r)

# 3.4 Fit Realized GARCH, supplying realizedVol
fit_real <- ugarchfit(
  spec         = spec_realgarch,
  data         = df_xts$r,
  realizedVol  = df_xts$x
)

# 4) Extract full‐sample conditional σ² ----------------------
# sigma(fit) returns the conditional standard deviation series
v_sg   <- sigma(fit_sg  )^2
v_vix  <- sigma(fit_vix )^2
v_epu  <- sigma(fit_epu )^2
v_real <- sigma(fit_real)^2

# realized variance series
rv_full <- df_xts$x

# 5) Merge into one xts and plot -----------------------------
all_xts <- merge(
  Realized      = rv_full,
  `GARCH(1,1)`  = v_sg,
  `GARCH-X (VIX)` = v_vix,
  `GARCH-X (EPU)` = v_epu,
  RealGARCH     = v_real
)

# autoplot will facet each series
autoplot(all_xts) +
  labs(
    title = "Full‐Sample Realized Variance & Model‐Implied σ²",
    x     = "Date",
    y     = "Variance",
    color = "Series"
  ) +
  theme_minimal()

# ----------------------------------------------------------
# 6) Compute full-sample MSE and QLIKE for each model
# ----------------------------------------------------------

# Define the loss functions
mse_func   <- function(f, r) mean((f - r)^2, na.rm=TRUE)
qlike_func <- function(f, r) mean(log(f) + r/f, na.rm=TRUE)

# realized variance full series (rv_full) and model‐implied σ² (v_sg, v_vix, v_epu, v_real)
# (these objects were created in steps 4 and 5 above)

# Align lengths (they all share the same index, df_xts)
# Compute metrics
results_metrics <- data.frame(
  Model = c("GARCH(1,1)", "GARCH-X (VIX)", "GARCH-X (EPU)", "RealGARCH"),
  MSE   = c(
    mse_func(v_sg,   rv_full),
    mse_func(v_vix,  rv_full),
    mse_func(v_epu,  rv_full),
    mse_func(v_real, rv_full)
  ),
  QLIKE = c(
    qlike_func(v_sg,   rv_full),
    qlike_func(v_vix,  rv_full),
    qlike_func(v_epu,  rv_full),
    qlike_func(v_real, rv_full)
  )
)

# Print the table of results
print(results_metrics)

















# ----------------------------------------------------------
# DM tests over the FULL SAMPLE
#   Models: GARCH11, GARCH-X(VIX), GARCH-X(EPU), RealGARCH
#   Losses: MSE (power=2) and QLIKE (power=1)
# ----------------------------------------------------------

# 0) Make sure you have these in your workspace:
#    - rv_full : xts of realized variance over the entire sample
#    - v_sg, v_vix, v_epu, v_real : xts of σ²_t from full-sample fits

# 1) Load forecast for dm.test()
library(forecast)

# 2) Extract coredata vectors (they all share the same dates/index)
real   <- coredata(rv_full)
sg     <- coredata(v_sg)
vix    <- coredata(v_vix)
epu    <- coredata(v_epu)
realg  <- coredata(v_real)

# 3) Put them in a named list for easy looping
forecasts_full <- list(
  GARCH11     = sg,
  GARCHX_VIX  = vix,
  GARCHX_EPU  = epu,
  RealGARCH   = realg
)
models <- names(forecasts_full)

# 4) Compute full‐sample loss series
loss_full_mse   <- lapply(forecasts_full, function(f) (f - real)^2)
loss_full_qlike <- lapply(forecasts_full, function(f) log(f) + real / f)

# 5) All unordered pairs of models
pairs <- combn(models, 2, simplify = FALSE)

# 6) Run DM tests pairwise
results_full <- data.frame(
  Model1  = character(),
  Model2  = character(),
  Metric  = character(),
  DM_stat = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(pr in pairs) {
  m1 <- pr[1]; m2 <- pr[2]
  
  # MSE-based DM test
  dm_mse <- dm.test(
    loss_full_mse[[m1]],
    loss_full_mse[[m2]],
    h     = 1,
    power = 2
  )
  
  # QLIKE-based DM test
  dm_ql <- dm.test(
    loss_full_qlike[[m1]],
    loss_full_qlike[[m2]],
    h     = 1,
    power = 1
  )
  
  # collect
  results_full <- rbind(
    results_full,
    data.frame(
      Model1  = m1, Model2 = m2,
      Metric  = "MSE",
      DM_stat = as.numeric(dm_mse$statistic),
      p_value = dm_mse$p.value,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Model1  = m1, Model2 = m2,
      Metric  = "QLIKE",
      DM_stat = as.numeric(dm_ql$statistic),
      p_value = dm_ql$p.value,
      stringsAsFactors = FALSE
    )
  )
}

# 7) Show results
print(results_full)




