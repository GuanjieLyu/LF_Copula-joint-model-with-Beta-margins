#################################
## This is the R script for copula-beta joint model in Lyu and Feng: Joint Copula-Based Beta 
## Regression with Covariate-Dependent Dispersion for Malaria Indicators
################################

library(MASS)     # mvrnorm
library(mgcv)     # gam, betar
library(GJRM)     # gjrm, copula.prob
library(dplyr)
library(purrr)
library(tidyr)
library(copula)
library(tibble)
library(ggplot2)
##========================================
## 1. Data-generating process 
##========================================
simulate_bivariate_beta_spatial <- function(
    n,
    tau_eval,
    beta_year  = c(-1.0, -0.8),
    gamma_disp = c(-1.6, -0.4, -1, -1.8),
    sd_spatial = 0.03,      # amplitude of spatial effect on eta-scale
    k_spatial  = 50,       # basis dimension (match your fitted k)
    coord_min  = 0,
    coord_max  = 10,
    seed       = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(requireNamespace("mgcv", quietly = TRUE))
  stopifnot(requireNamespace("copula", quietly = TRUE))
  stopifnot(requireNamespace("tibble", quietly = TRUE))
  
  ## 1) covariates
  x <- runif(n, -1, 1)
  z <- rbinom(n, 1, 0.5)
  
  ## 2) locations
  lat  <- runif(n, coord_min, coord_max)
  long <- runif(n, coord_min, coord_max)
  dat_loc <- data.frame(long = long, lat = lat)
  
  ## 3) build the SAME TP basis you will fit
  sm <- mgcv::s(long, lat, bs = "tp", k = k_spatial)
  sc <- mgcv::smoothCon(sm, data = dat_loc)[[1]]
  X  <- sc$X                 # n x M
  M  <- ncol(X)
  
  ## 4) simulate two independent spatial surfaces (NOT shared)
  beta1_sp <- rnorm(M)
  beta2_sp <- rnorm(M)
  b1_raw <- drop(X %*% beta1_sp)
  b2_raw <- drop(X %*% beta2_sp)
  
  ## scale to have sd = sd_spatial and mean 0 (eta-scale)
  b1 <- sd_spatial * as.numeric(scale(b1_raw))
  b2 <- sd_spatial * as.numeric(scale(b2_raw))
  
  ## 5) mean models + spatial
  f1 <- 1.2 * sin(pi * x); f1 <- f1 - mean(f1)
  f2 <- 0.8 * x^2;         f2 <- f2 - mean(f2)
  
  g  <- 0.03 * cos(pi * x); g <- g - mean(g)
  
  eta1 <- -0.7 + f1 + beta_year[1] * z + b1
  eta2 <- -0.5 + f2 + beta_year[2] * z + b2
  
  mu1 <- plogis(eta1)
  mu2 <- plogis(eta2)
  
  ## 6) dispersion (sigma with logit link)
  lp_sig1 <- gamma_disp[1] + gamma_disp[2] * z
  lp_sig2 <- gamma_disp[3] + g
  
  sigma1 <- plogis(lp_sig1)
  sigma2 <- plogis(lp_sig2)
  
  ## 7) sigma -> phi -> (alpha,beta)
  phi1 <- (1 - sigma1^2) / sigma1^2
  phi2 <- (1 - sigma2^2) / sigma2^2
  
  alpha1 <- mu1 * phi1
  beta1  <- (1 - mu1) * phi1
  
  alpha2 <- mu2 * phi2
  beta2  <- (1 - mu2) * phi2
  
  ## 8) copula uniforms
  theta_eval <- iTau(claytonCopula(), tau_eval)
  U <- copula::rCopula(n, copula::claytonCopula(param = theta_eval, dim = 2))
  
  ## 9) responses
  Y1 <- qbeta(U[,1], alpha1, beta1)
  Y2 <- qbeta(U[,2], alpha2, beta2)
  
  tibble::tibble(
    Y1 = Y1, Y2 = Y2,
    x = x, z = factor(z),
    long = long, lat = lat,
    b1 = b1, b2 = b2
  )
}


##==========================================================
## Model fitting functions
##==========================================================
fit_models <- function(dat) {
  
  # Joint copula – constant dispersion
  eq1 <- Y1 ~ s(x) + z + s(long, lat, bs = "tp")
  eq2 <- Y2 ~ s(x) + z + s(long, lat, bs = "tp")
  eq1.sig <- ~ 1
  eq2.sig <- ~ 1
  eq_theta <- ~ 1
  
  joint_const <- gjrm(
    formula = list(eq1, eq2, eq1.sig, eq2.sig, eq_theta),
    data = dat,
    margins = c("BE", "BE"),
    copula = "C0",
    model = "B"
  )
  
  # Independent GAMs (constant dispersion)
  gam1_const <- gamlss(list(eq1, eq1.sig), data = dat, family = "BE")
  gam2_const <- gamlss(list(eq2, eq2.sig), data = dat, family = "BE")
  
  # Joint copula – varying dispersion (TRUE model)
  eq1.sig <- ~ z
  eq2.sig <- ~ s(x)
  
  joint_var<- gjrm(
    formula = list(eq1, eq2, eq1.sig, eq2.sig, eq_theta),
    data = dat,
    margins = c("BE", "BE"),
    copula = "C0",
    model = "B"
  )
  
  gam1_var <- gamlss(list(eq1, eq1.sig), data = dat, family = "BE")
  gam2_var <- gamlss(list(eq2, eq2.sig), data = dat, family = "BE")
  
  list(
    gam1_const = gam1_const,
    gam2_const = gam2_const,
    joint_const = joint_const,
    gam1_var = gam1_var,
    gam2_var = gam2_var,
    joint_var   = joint_var
  )
}


true_upper_prob <- function(dat, c1, c2) {
  mean(dat$Y1 >= c1 & dat$Y2 >= c2)
}


independent_upper <- function(fit1, fit2, dat, c1, c2) {
  
  mu1 <- plogis(predict(fit1, dat, eq = 1, type = "link"))
  mu2 <- plogis(predict(fit2, dat, eq = 1, type = "link"))
  
  sigma1 <- plogis(predict(fit1, dat, eq = 2, type = "link"))
  sigma2 <- plogis(predict(fit2, dat, eq = 2, type = "link"))
  
  phi1 <- (1 - sigma1^2)/sigma1^2
  phi2 <- (1 - sigma2^2)/sigma2^2
  
  p1 <- 1 - pbeta(c1, mu1 * phi1, (1 - mu1) * phi1)
  p2 <- 1 - pbeta(c2, mu2 * phi2, (1 - mu2) * phi2)
  
  mean(p1 * p2)
}


# Copula-based upper-tail joint probability (simulation)
copula_upper <- function(fit, dat, c1, c2) {
  
  p1 <- sapply(seq_len(nrow(dat)), function(i) {
    copula.prob(
      x = fit, y1 = c1, y2 = 1,
      newdata = dat[i, ],
      joint = TRUE, cond = 0
    )$p12
  })
  
  p2 <- sapply(seq_len(nrow(dat)), function(i) {
    copula.prob(
      x = fit, y1 = 1, y2 = c2,
      newdata = dat[i, ],
      joint = TRUE, cond = 0
    )$p12
  })
  
  p12 <- sapply(seq_len(nrow(dat)), function(i) {
    copula.prob(
      x = fit, y1 = c1, y2 = c2,
      newdata = dat[i, ],
      joint = TRUE, cond = 0
    )$p12
  })
  
  mean(1 - p1 - p2 + p12)
}

##===================================================
## Run for replicates
##===================================================
set.seed(123)

n      <- 400    # sample size 
n_rep  <- 100    # number of replicates     
tau_evals <- c(0.1, 0.3, 0.6, 0.9)
results_list <- list()
counter <- 1

for (tau_eval in tau_evals) {
  for (rep in seq_len(n_rep)) {
    
    cat(sprintf("tau_eval = %.1f | rep %d / %d\n", tau_eval, rep, n_rep))
    
    dat <- simulate_bivariate_beta_spatial(n, tau_eval, seed = 2026 + rep + (tau_eval + 1) * 10)
    
    fits <- tryCatch(
      fit_models(dat),
      error = function(e) NULL
    )
    
    if (!is.null(fits)) {
      results_list[[counter]] <- list(
        tau_eval  = tau_eval,
        rep  = rep,
        dat  = dat,
        fits = fits
      )
      counter <- counter + 1
    }
  }
}

##=============================================
## filtering non-converged rep
##=============================================
conv_check_info <- function(fit) {
  # Capture printed output
  out <- capture.output(
    conv.check(fit)
  )
  
  # 1) Gradient line
  grad_line <- grep("Maximum absolute gradient value", out, value = TRUE)
  if (length(grad_line) == 0L) {
    grad <- NA_real_
  } else {
    # Extract numeric part after the colon
    grad_str <- sub(".*Maximum absolute gradient value:\\s*", "", grad_line[1])
    grad <- suppressWarnings(as.numeric(grad_str))
  }
  
  # 2) Positive definite flag
  pd <- any(grepl("Observed information matrix is positive definite", out))
  
  list(grad = grad, pd = pd)
}

info <- conv_check_info(results_list[[1]]$fits$joint_var)
info

conv_ok_one <- function(fit, grad_tol = 1e-1) {
  info <- conv_check_info(fit)
  (!is.na(info$grad)) && info$grad < grad_tol && isTRUE(info$pd)
}


is_converged <- function(res) {
  fits <- res$fits
  
  all(vapply(
    list(fits$gam1_const, fits$gam2_const,
         fits$gam1_var, fits$gam2_var,
         fits$joint_const, fits$joint_var),
    conv_ok_one,
    logical(1)
  ))
}

is_ok <- vapply(results_list, is_converged, logical(1))
table(is_ok) 
results_ok <- results_list[is_ok]

which(!is_ok)
conv.check(results_list[[8]]$fits$joint_var)


#######################################
# Safe AIC/BIC
#######################################
safe_AIC <- function(x) {
  if (is.null(x)) return(NA_real_)
  out <- tryCatch(AIC(x), error = function(e) NA_real_)
  as.numeric(out)
}

safe_BIC <- function(x) {
  if (is.null(x)) return(NA_real_)
  out <- tryCatch(BIC(x), error = function(e) NA_real_)
  as.numeric(out)
}

#######################################
# Compute AIC/BIC for each replicate
# (robust to list-of-results or list-of-lists)
#######################################
n <- length(results_ok)

ic_df <- data.frame(
  tau_eval         = numeric(n),
  rep              = integer(n),
  AIC_indep_const  = numeric(n),
  AIC_joint_const  = numeric(n),
  AIC_indep_var    = numeric(n),
  AIC_joint_var    = numeric(n),
  BIC_indep_const  = numeric(n),
  BIC_joint_const  = numeric(n),
  BIC_indep_var    = numeric(n),
  BIC_joint_var    = numeric(n)
)

for (i in seq_len(n)) {
  x <- results_ok[[i]]
  fits <- x$fits
  
  ic_df$tau_eval[i] <- x$tau_eval
  ic_df$rep[i]      <- x$rep
  
  ic_df$AIC_indep_const[i] <-
    safe_AIC(fits$gam1_const) + safe_AIC(fits$gam2_const)
  
  ic_df$AIC_joint_const[i] <-
    safe_AIC(fits$joint_const)
  
  ic_df$AIC_indep_var[i] <-
    safe_AIC(fits$gam1_var) + safe_AIC(fits$gam2_var)
  
  ic_df$AIC_joint_var[i] <-
    safe_AIC(fits$joint_var)
  
  ic_df$BIC_indep_const[i] <-
    safe_BIC(fits$gam1_const) + safe_BIC(fits$gam2_const)
  
  ic_df$BIC_joint_const[i] <-
    safe_BIC(fits$joint_const)
  
  ic_df$BIC_indep_var[i] <-
    safe_BIC(fits$gam1_var) + safe_BIC(fits$gam2_var)
  
  ic_df$BIC_joint_var[i] <-
    safe_BIC(fits$joint_var)
}

#######################################
# Average AIC/BIC over replicates
#######################################
ic_summary <- ic_df %>%
  group_by(tau_eval) %>%
  summarise(
    across(starts_with("AIC"), mean, na.rm = TRUE),
    across(starts_with("BIC"), mean, na.rm = TRUE),
    .groups = "drop"
  )

#######################################
# LaTeX table (xtable)
#######################################
library(xtable)

xtable(
  ic_summary,
  digits  = c(0, 1, rep(3, ncol(ic_summary) - 1)),
  caption = "Average AIC and BIC across simulation replicates",
  label   = "tab:sim_aicbic"
)


##############################
# Estimated tau_values
#############################
tau_param_df <- map_dfr(results_ok, function(x) {
  
  fit_c <- x$fits$joint_const
  fit_v <- x$fits$joint_var
  
  tibble(
    n       = if (!is.null(x$n)) x$n else nrow(x$dat),  # keep if available
    tau_eval = x$tau_eval,
    rep      = x$rep,
    
    tau_hat_jc = if (!is.null(fit_c)) tau(claytonCopula(fit_c$theta[1])) else NA_real_,
    tau_hat_jv = if (!is.null(fit_v)) tau(claytonCopula(fit_v$theta[1])) else NA_real_
  )
})

tau_summary <- tau_param_df %>%
  pivot_longer(
    cols = c(tau_hat_jc, tau_hat_jv),
    names_to = "model",
    values_to = "tau_hat"
  ) %>%
  mutate(
    model = recode(model,
                   tau_hat_jc = "Joint constant",
                   tau_hat_jv = "Joint varying"
    )
  ) %>%
  group_by(tau_eval, model) %>%
  summarise(
    bias = mean(tau_hat - tau_eval, na.rm = TRUE),
    rmse = sqrt(mean((tau_hat - tau_eval)^2, na.rm = TRUE)),
    n_rep = sum(!is.na(tau_hat)),
    .groups = "drop"
  )

tau_summary <- tau_summary %>%
  mutate(across(c(bias, rmse), ~ round(.x, 3)))

######################################
# z-coefficients
######################################
extract_z_gamlss <- function(fit, margin = c("mu", "sigma")) {
  margin <- match.arg(margin)
  if (is.null(fit)) return(NA_real_)
  
  s <- summary(fit)
  
  tab <- switch(
    margin,
    mu = s$tableP1,
    sigma = s$tableP2,
  )
  
  # row name is "z1" because z is a factor
  if (!"z1" %in% rownames(tab)) return(NA_real_)
  
  tab["z1", "Estimate"]
}

extract_z_gjrm <- function(fit, margin = c("mu1", "mu2", "sigma1")) {
  margin <- match.arg(margin)
  if (is.null(fit)) return(NA_real_)
  
  s <- summary(fit)
  
  tab <- switch(
    margin,
    mu1 = s$tableP1,
    mu2 = s$tableP2,
    sigma1 = s$tableP3
  )
  
  # row name is "z1" because z is a factor
  if (!"z1" %in% rownames(tab)) return(NA_real_)
  
  tab["z1", "Estimate"]
}



library(purrr)
library(dplyr)
library(tibble)

z_param_df <- map_dfr(results_ok, function(x) {
  
  fits <- x$fits
  
  tibble(
    tau_eval = x$tau_eval,
    rep = x$rep,
    
    # separate constant models
    gam_const_mu1_z = extract_z_gamlss(fits$gam1_const, "mu"),
    gam_const_mu2_z = extract_z_gamlss(fits$gam2_const, "mu"),
    
    # separate varying models
    gam_mu1_z = extract_z_gamlss(fits$gam1_var, "mu"),
    gam_sigma1_z = extract_z_gamlss(fits$gam1_var, "sigma"),
    gam_mu2_z = extract_z_gamlss(fits$gam2_var, "mu"),
    
    # Joint copula – constant dispersion
    jc_mu1_z  = extract_z_gjrm(fits$joint_const, "mu1"),
    jc_mu2_z  = extract_z_gjrm(fits$joint_const, "mu2"),
    
    # Joint copula – varying dispersion
    jv_mu1_z  = extract_z_gjrm(fits$joint_var, "mu1"),
    jv_sigma1_z  = extract_z_gjrm(fits$joint_var, "sigma1"),
    jv_mu2_z  = extract_z_gjrm(fits$joint_var, "mu2")
  )
})


true_mu1_z <- c(-1.0, -0.8)[1]
true_mu2_z <- c(-1.0, -0.8)[2]
true_sigma1_z <- -0.4

bias_rmse_z <- z_param_df %>%
  group_by(tau_eval) %>%
  summarise(
    
    ## ---- Y1 ----
    bias_gam_const_mu1 = mean(gam_const_mu1_z - true_mu1_z, na.rm = TRUE),
    rmse_gam_const_mu1 = sqrt(mean((gam_const_mu1_z - true_mu1_z)^2, na.rm = TRUE)),
    
    bias_gam_mu1 = mean(gam_mu1_z - true_mu1_z, na.rm = TRUE),
    rmse_gam_mu1 = sqrt(mean((gam_mu1_z - true_mu1_z)^2, na.rm = TRUE)),
    
    bias_gam_sigma1 = mean(gam_sigma1_z - true_sigma1_z, na.rm = TRUE),
    rmse_gam_sigma1 = sqrt(mean((gam_sigma1_z - true_sigma1_z)^2, na.rm = TRUE)),
    
    
    bias_jc_mu1  = mean(jc_mu1_z - true_mu1_z, na.rm = TRUE),
    rmse_jc_mu1  = sqrt(mean((jc_mu1_z - true_mu1_z)^2, na.rm = TRUE)),
    
    bias_jv_mu1  = mean(jv_mu1_z - true_mu1_z, na.rm = TRUE),
    rmse_jv_mu1  = sqrt(mean((jv_mu1_z - true_mu1_z)^2, na.rm = TRUE)),
    
    bias_jv_sigma1  = mean(jv_sigma1_z - true_sigma1_z, na.rm = TRUE),
    rmse_jv_sigma1  = sqrt(mean((jv_sigma1_z - true_sigma1_z)^2, na.rm = TRUE)),
    
    ## ---- Y2 ----
    bias_gam_const_mu2 = mean(gam_const_mu2_z - true_mu2_z, na.rm = TRUE),
    rmse_gam_const_mu2 = sqrt(mean((gam_const_mu2_z - true_mu2_z)^2, na.rm = TRUE)),
    
    bias_gam_mu2 = mean(gam_mu2_z - true_mu2_z, na.rm = TRUE),
    rmse_gam_mu2 = sqrt(mean((gam_mu2_z - true_mu2_z)^2, na.rm = TRUE)),
    
    bias_jc_mu2  = mean(jc_mu2_z - true_mu2_z, na.rm = TRUE),
    rmse_jc_mu2  = sqrt(mean((jc_mu2_z - true_mu2_z)^2, na.rm = TRUE)),
    
    bias_jv_mu2  = mean(jv_mu2_z - true_mu2_z, na.rm = TRUE),
    rmse_jv_mu2  = sqrt(mean((jv_mu2_z - true_mu2_z)^2, na.rm = TRUE)),
    
    .groups = "drop"
  )

bias_rmse_mu1 <- bias_rmse_z %>%
  dplyr::select(
    tau_eval,
    bias_gam_const_mu1, bias_gam_mu1, bias_jc_mu1, bias_jv_mu1,
    rmse_gam_const_mu1, rmse_gam_mu1, rmse_jc_mu1, rmse_jv_mu1
  ) %>%
  rename(
    Bias_GAM_CONST  = bias_gam_const_mu1,
    Bias_GAM  = bias_gam_mu1,
    Bias_JointConst = bias_jc_mu1,
    Bias_JointVar   = bias_jv_mu1,
    RMSE_GAM_CONST  = rmse_gam_const_mu1,
    RMSE_GAM  = rmse_gam_mu1,
    RMSE_JointConst = rmse_jc_mu1,
    RMSE_JointVar   = rmse_jv_mu1
  )

bias_rmse_mu1 <- bias_rmse_mu1 %>%
  mutate(across(-tau_eval, ~ round(.x, 3)))
print(bias_rmse_mu1)
bias_rmse_mu1[, 8:9]
##-----------------------------------------------------

bias_rmse_mu2 <- bias_rmse_z %>%
  dplyr::select(
    tau_eval,
    bias_gam_const_mu2, bias_gam_mu2, bias_jc_mu2, bias_jv_mu2,
    rmse_gam_const_mu2, rmse_gam_mu2, rmse_jc_mu2, rmse_jv_mu2
  ) %>%
  rename(
    Bias_GAM_CONST  = bias_gam_const_mu2,
    Bias_GAM  = bias_gam_mu2,
    Bias_JointConst = bias_jc_mu2,
    Bias_JointVar   = bias_jv_mu2,
    RMSE_GAM_CONST  = rmse_gam_const_mu2,
    RMSE_GAM  = rmse_gam_mu2,
    RMSE_JointConst = rmse_jc_mu2,
    RMSE_JointVar   = rmse_jv_mu2
  )

bias_rmse_mu2 <- bias_rmse_mu2 %>%
  mutate(across(-tau_eval, ~ round(.x, 3)))
print(bias_rmse_mu2)
bias_rmse_mu2[, 8:9]
##-------------------------------------------------------------
bias_rmse_sigma1 <- bias_rmse_z %>%
  dplyr::select(
    tau_eval,
    bias_gam_sigma1, bias_jv_sigma1,
    rmse_gam_sigma1, rmse_jv_sigma1
  ) %>%
  rename(
    Bias_GAM  = bias_gam_sigma1,
    Bias_JointVar   = bias_jv_sigma1,
    RMSE_GAM  = rmse_gam_sigma1,
    RMSE_JointVar   = rmse_jv_sigma1
  )

bias_rmse_sigma1<- bias_rmse_sigma1 %>%
  mutate(across(-tau_eval, ~ round(.x, 3)))
print(bias_rmse_sigma1)

##------------------------------------------------------------
library(xtable)

xt <- xtable(
  bias_rmse_mu1,
  caption = "Bias and RMSE of $z$ coefficient in the mean model for $Y_1$",
  label = "tab:bias_rmse_mu1",
  digits = c(0, 1, rep(3, ncol(bias_rmse_mu1) - 1))
)

print(
  xt,
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity
)




####################################
## box_plots
####################################
z_param_df400 <- readRDS("z_param_df400.rds")
z_param_df800 <- readRDS("z_param_df800.rds")
z_param_df1600 <- readRDS("z_param_df1600.rds")

str(z_param_df400)
str(z_param_df800)
str(z_param_df1600)

##===============================================================
df400  <- z_param_df400  %>% mutate(SampleSize = 400)
df800  <- z_param_df800  %>% mutate(SampleSize = 800)
df1600 <- z_param_df1600 %>% mutate(SampleSize = 1600)

sim_data <- bind_rows(df400, df800, df1600)

sim_mu1 <- sim_data %>%
  select(
    tau_eval, SampleSize,
    gam_const_mu1_z,
    gam_mu1_z,
    jc_mu1_z,
    jv_mu1_z
  ) %>%
  pivot_longer(
    cols = c(gam_const_mu1_z, gam_mu1_z, jc_mu1_z, jv_mu1_z),
    names_to  = "Model",
    values_to = "BetaHat"
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = c("gam_const_mu1_z", "gam_mu1_z",
                 "jc_mu1_z", "jv_mu1_z"),
      labels = c("IN-CD", "IN-VD",
                 "JC-CD", "JC-VD")
    ),
    SampleSize = factor(SampleSize),
    Tau = tau_eval
  )

ggplot(sim_mu1, aes(x = SampleSize, y = BetaHat, fill = Model)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65,
               position = position_dodge(width = 0.75)) +
  facet_wrap(~ Tau, nrow = 2,
             labeller = label_bquote(tau == .(Tau))) +
  scale_fill_manual(values = c(
    "IN-CD" = "#8DB4E2",
    "IN-VD"  = "#F4C27A",
    "JC-CD"  = "#A8D5A2",
    "JC-VD"   = "#E7A1B0"
  )) +
  geom_hline(yintercept = -1.0, linetype = "dashed", linewidth = 0.4) +
  theme_bw(base_size = 13) +
  labs(
    title = bquote("Estimated " * beta[1,z] * " (Outcome 1 mean component)"),
    x = "Sample size",
    y = expression(hat(beta)[1,z]),
    fill = "Model"
  ) + 
  theme(
    legend.position = "bottom")
##===========================================================================

sim_mu2 <- sim_data %>%
  select(
    tau_eval, SampleSize,
    gam_const_mu2_z,
    gam_mu2_z,
    jc_mu2_z,
    jv_mu2_z
  ) %>%
  pivot_longer(
    cols = c(gam_const_mu2_z, gam_mu2_z, jc_mu2_z, jv_mu2_z),
    names_to  = "Model",
    values_to = "BetaHat"
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = c("gam_const_mu2_z", "gam_mu2_z",
                 "jc_mu2_z", "jv_mu2_z"),
      labels = c("IN-CD", "IN-VD",
                 "JC-CD", "JC-VD")
    ),
    SampleSize = factor(SampleSize),
    Tau = tau_eval
  )

ggplot(sim_mu2, aes(x = SampleSize, y = BetaHat, fill = Model)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65,
               position = position_dodge(width = 0.75)) +
  facet_wrap(~ Tau, nrow = 2,
             labeller = label_bquote(tau == .(Tau))) +
  scale_fill_manual(values = c(
    "IN-CD" = "#8DB4E2",
    "IN-VD"  = "#F4C27A",
    "JC-CD"  = "#A8D5A2",
    "JC-VD"   = "#E7A1B0"
  )) +
  geom_hline(yintercept = -0.8, linetype = "dashed", linewidth = 0.4) +
  theme_bw(base_size = 13) +
  labs(
    title = bquote("Estimated " * beta[2,z] * " (Outcome 2 mean component)"),
    x = "Sample size",
    y = expression(hat(beta)[2,z]),
    fill = "Model"
  ) + 
  theme(
    legend.position = "bottom")
##=============================================================

sim_sigma1 <- sim_data %>%
  select(
    tau_eval, SampleSize,
    gam_sigma1_z,
    jv_sigma1_z
  ) %>%
  pivot_longer(
    cols = c(gam_sigma1_z, jv_sigma1_z),
    names_to  = "Model",
    values_to = "BetaHat"
  ) %>%
  mutate(
    Model = factor(
      Model,
      levels = c("gam_sigma1_z", "jv_sigma1_z"),
      labels = c("IN-VD", "JC-VD")
    ),
    SampleSize = factor(SampleSize),
    Tau = tau_eval
  )

ggplot(sim_sigma1, aes(x = SampleSize, y = BetaHat, fill = Model)) +
  geom_boxplot(outlier.size = 0.7, width = 0.65,
               position = position_dodge(width = 0.75)) +
  facet_wrap(~ Tau, nrow = 2,
             labeller = label_bquote(tau == .(Tau))) +
  scale_fill_manual(values = c(
    "IN-VD" = "#F4C27A",
    "JC-VD" = "#E7A1B0"
  )) +
  geom_hline(yintercept = -0.4, linetype = "dashed", linewidth = 0.4) +
  theme_bw(base_size = 13) +
  labs(
    title = bquote("Estimated " * gamma[1,z] * " (Outcome 1 dispersion component)"),
    x = "Sample size",
    y = expression(hat(gamma)[1,z]),
    fill = "Model"
  ) + 
  theme(
    legend.position = "bottom")
##===========================================================
## Boxplots in one row
##===========================================================
library(dplyr)
library(tidyr)
library(ggplot2)

tau_levels <- c(0.1, 0.3, 0.6, 0.9)

model_levels4 <- c("gam_const", "gam_var", "jc_const", "jc_var")

df_mu1 <- sim_data %>%
  transmute(
    Tau = tau_eval,
    SampleSize = factor(SampleSize),
    Parameter = "beta1z",
    gam_const = gam_const_mu1_z,
    gam_var   = gam_mu1_z,
    jc_const  = jc_mu1_z,
    jc_var    = jv_mu1_z
  ) %>%
  pivot_longer(all_of(model_levels4), names_to = "Model", values_to = "Estimate")

df_mu2 <- sim_data %>%
  transmute(
    Tau = tau_eval,
    SampleSize = factor(SampleSize),
    Parameter = "beta2z",
    gam_const = gam_const_mu2_z,
    gam_var   = gam_mu2_z,
    jc_const  = jc_mu2_z,
    jc_var    = jv_mu2_z
  ) %>%
  pivot_longer(all_of(model_levels4), names_to = "Model", values_to = "Estimate")

df_sig1 <- sim_data %>%
  transmute(
    Tau = tau_eval,
    SampleSize = factor(SampleSize),
    Parameter = "gamma1z",
    gam_var = gam_sigma1_z,
    jc_var  = jv_sigma1_z
  ) %>%
  pivot_longer(c(gam_var, jc_var), names_to = "Model", values_to = "Estimate")

plot_df <- bind_rows(df_mu1, df_mu2, df_sig1) %>%
  mutate(
    Tau = factor(Tau, levels = tau_levels),
    Model = factor(
      Model,
      levels = c("gam_const", "gam_var", "jc_const", "jc_var"),
      labels = c("IN-CD", "IN-VD", "JC-CD", "JC-VD")
    ),
    # Use parseable math strings directly as factor labels
    Parameter = factor(
      Parameter,
      levels = c("beta1z", "beta2z", "gamma1z"),
      labels = c("hat(beta)[1,z]", "hat(beta)[2,z]", "hat(gamma)[1,z]")
    )
  )


# Define true values for each parameter (adjust these to your actual true values)
true_values <- data.frame(
  Parameter = c("param1", "param2", "param3"),  # match your Parameter levels exactly
  TrueValue = c(-1.0, -0.8, -0.4)                 # replace with your actual true values
)

ggplot(plot_df, aes(x = SampleSize, y = Estimate, fill = Model)) +
  geom_boxplot(
    outlier.size = 0.7,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +
  facet_grid(
    Parameter ~ Tau,
    scales = "free_y",
    labeller = labeller(
      Tau = function(values) paste0("\u03C4 = ", values),
      Parameter = label_parsed
    )
  ) +
  scale_fill_manual(values = c(
    "IN-CD" = "#8DB4E2",
    "IN-VD" = "#F4C27A",
    "JC-CD" = "#A8D5A2",
    "JC-VD" = "#E7A1B0"
  ), drop = FALSE) +
  geom_hline(
    data = data.frame(Parameter = unique(plot_df$Parameter)[1], yintercept = -1.0),
    aes(yintercept = yintercept),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    data = data.frame(Parameter = unique(plot_df$Parameter)[2], yintercept = -0.8),
    aes(yintercept = yintercept),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    data = data.frame(Parameter = unique(plot_df$Parameter)[3], yintercept = -0.4),
    aes(yintercept = yintercept),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "Sample size",
    y = NULL,
    fill = "Model"
  ) +  theme(
    legend.position = "bottom",
    panel.spacing.x = unit(0.9, "lines"),
    panel.spacing.y = unit(0.9, "lines"),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12, angle = 0),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    strip.background.x = element_rect(fill = "lightgrey", color = "black"),  # grey top labels
    strip.background.y = element_rect(fill = "lightgrey", color = "black")   # grey right labels
  )
  


