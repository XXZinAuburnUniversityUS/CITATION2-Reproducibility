##  demonstration code for censor–weight + standardized survival

library(dplyr)
library(survival)
library(splines)
library(boot)

## 1. Basic preprocessing: define chemo_binary and per-protocol censoring
prep_pp_data <- function(data) {
  data <- data %>%
    arrange(PATIENT_ID, cycle) %>%
    group_by(PATIENT_ID) %>%
    mutate(
      death_next = ifelse(is.na(lead(death)), 0, lead(death)),
      time       = row_number() - 1,
      chemo_binary = case_when(
        time == 0 & Treatment == "ICI_Chemo" ~ 1,
        time == 0 & Treatment == "ICI_Mono"  ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    mutate(
      chemo_binary = ifelse(any(chemo_binary == 1, na.rm = TRUE), 1,
                            ifelse(any(chemo_binary == 0, na.rm = TRUE), 0, NA_real_))
    ) %>%
    ungroup()
  
  mono  <- data %>% filter(chemo_binary == 0)
  chemo <- data %>% filter(chemo_binary == 1)
  
  mono_pp <- mono %>%
    group_by(PATIENT_ID) %>%
    mutate(
      first_non_mono = min(ifelse(cycle %in% 1:4 & Treatment != "ICI_Mono", cycle, Inf)),
      censor         = ifelse(cycle == first_non_mono & first_non_mono < Inf, 1L, 0L)
    ) %>%
    filter(cycle <= first_non_mono | first_non_mono == Inf) %>%
    ungroup() %>%
    mutate(arm = "ICI-Mono")
  
  chemo_pp <- chemo %>%
    group_by(PATIENT_ID) %>%
    mutate(
      first_non_chemo = min(ifelse(cycle %in% 1:4 & Treatment != "ICI_Chemo", cycle, Inf)),
      censor          = ifelse(cycle == first_non_chemo & first_non_chemo < Inf, 1L, 0L)
    ) %>%
    filter(cycle <= first_non_chemo | first_non_chemo == Inf) %>%
    ungroup() %>%
    mutate(arm = "ICI-Chemo")
  
  bind_rows(chemo_pp, mono_pp)
}

## 2. Fit censoring weights and discrete-time model, return main objects
run_ccw_demo <- function(data, weight_cycles = 2:4) {
  pp_dat <- prep_pp_data(data)
  
  # split by arm
  chemo_pp <- pp_dat %>% filter(chemo_binary == 1)
  mono_pp  <- pp_dat %>% filter(chemo_binary == 0)
  
  # simple covariate set for demonstration
  num_form <- (censor == 0) ~ ns(cycle, 2) + age + SEX_CATEGORY
  den_form <- (censor == 0) ~ ns(cycle, 2) + age + SEX_CATEGORY + Charlson
  
  # arm-specific weight models
  fit_weights <- function(dat) {
    model_dat <- dat %>% filter(cycle %in% weight_cycles)
    num <- glm(num_form, family = quasibinomial(), data = model_dat)
    den <- glm(den_form, family = quasibinomial(), data = model_dat)
    
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
  
  chemo_w <- fit_weights(chemo_pp) %>% mutate(chemo_binary = 1)
  mono_w  <- fit_weights(mono_pp)  %>% mutate(chemo_binary = 0)
  
  combined <- bind_rows(chemo_w, mono_w)
  
  ## 3. Weighted Cox model (clone–censor–weight)
  cox_fit <- coxph(
    Surv(cycle - 1, cycle, death_next) ~ chemo_binary + age + SEX_CATEGORY + Charlson,
    data    = combined,
    weights = sw_cum,
    robust  = TRUE,
    cluster = PATIENT_ID
  )
  
  ## 4. Weighted discrete-time hazard model
  disc_fit <- glm(
    death_next ~ ns(cycle, 4) * chemo_binary + age + SEX_CATEGORY + Charlson,
    family  = quasibinomial(),
    data    = combined,
    weights = sw_cum
  )
  
  ## 5. Standardized survival curves under each strategy (no bootstrap shown)
  base_dat <- combined %>%
    filter(cycle == 1) %>%
    distinct(PATIENT_ID, age, SEX_CATEGORY, Charlson)
  
  make_surv <- function(chemo_value, t_max = 48) {
    t_vec <- 1:t_max
    newdat <- base_dat[rep(seq_len(nrow(base_dat)), each = t_max), ]
    newdat$cycle       <- rep(t_vec, times = nrow(base_dat))
    newdat$chemo_binary <- chemo_value
    
    newdat$pred_h <- predict(disc_fit, newdata = newdat, type = "response")
    
    newdat <- newdat %>%
      group_by(PATIENT_ID) %>%
      arrange(cycle) %>%
      mutate(S = cumprod(1 - pred_h)) %>%
      ungroup()
    
    newdat %>%
      group_by(cycle) %>%
      summarise(S = mean(S), .groups = "drop")
  }
  
  surv_mono  <- make_surv(chemo_value = 0)
  surv_chemo <- make_surv(chemo_value = 1)
  
  list(
    data_pp      = combined,
    cox_summary  = summary(cox_fit),
    disc_summary = summary(disc_fit),
    surv_mono    = surv_mono,
    surv_chemo   = surv_chemo
  )
}
