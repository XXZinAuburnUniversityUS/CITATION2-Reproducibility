## -------------------------------------------------
## Assumes processed_data3 already exists:
## - long format, 1 row per patient-cycle
## - has: PATIENT_ID, cycle, death1, time, Treatment, chemo_binary, 
##         ICI, DIAGNOSIS_YEAR1, age, TTI, region, METRO_CATEGORY,
##         SEX_CATEGORY, SIMPLE_RACE_CATEGORY1, YOST, histology,
##         baseline_* covariates, Charlson, adherence_ED, etc.
## -------------------------------------------------

library(dplyr)
library(survival)
library(splines)
library(boot)

## 1. Define chemo vs mono cohorts and per-protocol (no-grace) censoring

# Identify baseline chemo patients
patient_chemo <- processed_data3 %>% 
  filter(cycle == 1, chemo_binary == 1) %>% 
  select(PATIENT_ID) %>% 
  distinct()

CHEMO <- processed_data3 %>%
  filter(PATIENT_ID %in% patient_chemo$PATIENT_ID)

MONO <- processed_data3 %>%
  filter(!PATIENT_ID %in% patient_chemo$PATIENT_ID)

# Per-protocol censoring for ICI-Mono arm (no grace)
MONO_pp <- MONO %>%
  group_by(PATIENT_ID) %>%
  mutate(
    first_non_mono_cycle = suppressWarnings(
      min(ifelse(cycle %in% 1:4 & Treatment != "ICI_Mono", cycle, Inf))
    ),
    has_deviation_1_4 = first_non_mono_cycle != Inf,
    censor = if_else(has_deviation_1_4 & cycle == first_non_mono_cycle, 1L, 0L)
  ) %>%
  filter(cycle <= first_non_mono_cycle | !has_deviation_1_4) %>%
  ungroup() %>%
  mutate(arm_name = "ICI Mono")

# Per-protocol censoring for ICI-Chemo arm (no grace)
CHEMO_pp <- CHEMO %>%
  group_by(PATIENT_ID) %>%
  mutate(
    first_non_chemo_cycle = suppressWarnings(
      min(ifelse(cycle %in% 1:4 & Treatment != "ICI_Chemo", cycle, Inf))
    ),
    has_deviation_1_4 = first_non_chemo_cycle != Inf,
    censor = if_else(has_deviation_1_4 & cycle == first_non_chemo_cycle, 1L, 0L)
  ) %>%
  filter(cycle <= first_non_chemo_cycle | !has_deviation_1_4) %>%
  ungroup() %>%
  mutate(arm_name = "ICI + Chemo")

CHEMO_pp_nograce <- CHEMO_pp
MONO_pp_nograce  <- MONO_pp

## -------------------------------------------------
## 2. Main demonstration function:
##    - CCW weighting in each arm
##    - weighted discrete-time model
##    - bootstrap standardized survival (chemo vs mono)
##    All formulas are passed in as arguments.
## -------------------------------------------------

run_ccw_demo <- function(
    data_chemo,       # per-protocol ICI-Chemo long data
    data_mono,        # per-protocol ICI-Mono long data
    censor_var,       # name of censor variable, e.g. "censor"
    weight_cycles,    # e.g., 2:4
    num_formula,      # numerator model for censoring weights
    den_formula,      # denominator model for censoring weights
    outcome_formula,  # discrete-time outcome model formula
    t_max   = 48,     # maximum cycle for standardized survival
    R_boot  = 300     # number of bootstrap replicates
) {
  # Helper to fit weights within one arm
  fit_weights <- function(dat) {
    model_dat <- dat %>% filter(cycle %in% weight_cycles)
    
    num <- glm(num_formula, family = quasibinomial(), data = model_dat)
    den <- glm(den_formula, family = quasibinomial(), data = model_dat)
    
    dat <- dat %>%
      mutate(
        p_num = ifelse(
          cycle %in% weight_cycles,
          predict(num, newdata = ., type = "response"),
          1
        ),
        p_den = ifelse(
          cycle %in% weight_cycles,
          predict(den, newdata = ., type = "response"),
          1
        ),
        sw_raw = p_num / p_den
      )
    
    trunc_lo <- quantile(dat$sw_raw, 0.025, na.rm = TRUE)
    trunc_hi <- quantile(dat$sw_raw, 0.975, na.rm = TRUE)
    
    dat %>%
      mutate(
        sw_trunc = pmin(pmax(sw_raw, trunc_lo), trunc_hi)
      ) %>%
      group_by(PATIENT_ID) %>%
      arrange(cycle) %>%
      mutate(sw_cum = cumprod(sw_trunc)) %>%
      ungroup()
  }
  
  # Fit weights in each arm separately
  chemo_w <- fit_weights(data_chemo) %>%
    mutate(chemo_binary = 1)
  mono_w  <- fit_weights(data_mono) %>%
    mutate(chemo_binary = 0)
  
  # Combine arms
  combined_pp <- bind_rows(chemo_w, mono_w)
  
  # Weighted discrete-time hazard model
  disc_fit <- glm(
    outcome_formula,
    family  = quasibinomial(),
    data    = combined_pp,
    weights = sw_cum
  )
  
  # Baseline covariates for standardization (1 row per patient at cycle 1)
  baseline <- combined_pp %>%
    filter(cycle == 1) %>%
    distinct(PATIENT_ID, .keep_all = TRUE)
  
  # Single-bootstrap standardized survival (returns c(S_mono, S_chemo))
  standardize_once <- function(indices) {
    base_b <- baseline[indices, ]
    
    expand_for_arm <- function(chemo_value) {
      t_vec <- 1:t_max
      newdat <- base_b[rep(seq_len(nrow(base_b)), each = t_max), ]
      newdat$cycle        <- rep(t_vec, times = nrow(base_b))
      newdat$chemo_binary <- chemo_value
      newdat$pred_h       <- predict(disc_fit, newdata = newdat, type = "response")
      
      newdat %>%
        group_by(PATIENT_ID) %>%
        arrange(cycle) %>%
        mutate(S = cumprod(1 - pred_h)) %>%
        ungroup() %>%
        group_by(cycle) %>%
        summarise(S = mean(S), .groups = "drop") %>%
        pull(S)
    }
    
    S_mono  <- expand_for_arm(chemo_value = 0)
    S_chemo <- expand_for_arm(chemo_value = 1)
    
    c(S_mono, S_chemo)
  }
  
  # Bootstrap over individuals
  set.seed(123)
  boot_res <- boot(
    data      = baseline,
    statistic = function(d, i) standardize_once(i),
    R         = R_boot
  )
  
  # Pack point estimates at each cycle (t0)
  num_t <- t_max
  S_mono_t0  <- boot_res$t0[1:num_t]
  S_chemo_t0 <- boot_res$t0[(num_t + 1):(2 * num_t)]
  
  surv_estimates <- data.frame(
    cycle     = 1:num_t,
    S_mono    = S_mono_t0,
    S_chemo   = S_chemo_t0,
    diff_S    = S_chemo_t0 - S_mono_t0
  )
  
  invisible(list(
    combined_pp      = combined_pp,
    outcome_model    = summary(disc_fit),
    boot_survival    = boot_res,
    surv_estimates   = surv_estimates
  ))
}

## -------------------------------------------------
## 3. Example usage (formulas intentionally left for users to specify)
##    In the supplement / GitHub, you can show them as placeholders.
## -------------------------------------------------

# num_formula     <- (censor == 0) ~ ns(cycle, 2) + [user-specified covariates]
# den_formula     <- (censor == 0) ~ ns(cycle, 2) + [user-specified covariates + additional predictors]
# outcome_formula <- death1        ~ ns(cycle, 4) * chemo_binary + [user-specified covariates]

# res_demo <- run_ccw_demo(
#   data_chemo      = CHEMO_pp_nograce,
#   data_mono       = MONO_pp_nograce,
#   censor_var      = "censor",
#   weight_cycles   = 2:4,
#   num_formula     = num_formula,
#   den_formula     = den_formula,
#   outcome_formula = outcome_formula,
#   t_max           = 48,
#   R_boot          = 300
# )
