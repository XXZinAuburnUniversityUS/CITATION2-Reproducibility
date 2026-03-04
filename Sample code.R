###############################################################################
# CITATION-2: Per-Protocol Analysis with One-Cycle Grace Period
# - Data: SEER-Medicare processed longitudinal dataset (cycle-level; 1 cycle = 21 days)
# - Objective: Compare continuation strategy using IPW + discrete-time hazard
# - Outputs:
#   (1) Censoring flow plot
#   (2) Stabilized weight distributions
#   (3) Cox models (0–144, 0–24, 0–48 weeks)
#   (4) Discrete-time hazard model
#   (5) Bootstrap-standardized survival curves + risk difference table (DOCX)
#   (6) Model summaries (DOCX)
###############################################################################

## ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(glue)
  library(splines)
  library(boot)
  library(survival)
  library(scales)
  library(flextable)
  library(officer)
  library(patchwork)
})

## ---- Housekeeping ----
rm(list = ls())

## ---- I/O Paths ----
data_path <- "M:/Xiangzhong Xue [XZX0022]/seer_medicare practice/r code/rdata/processed_data_CCI.rds"

save_path <- "M:/Xiangzhong Xue [XZX0022]/seer_medicare practice/r code/Aim 1_time varing/Revision/Revision 2"

## ---- Load Data ----
processed_data_raw <- readRDS(data_path)

## ---- Select Variables ----
vars_to_keep <- c(
  # Model variables
  "death", "cycle", "ICI", "DIAGNOSIS_YEAR1", "age", "TTI",
  "REGION", "METRO_CATEGORY", "SEX_CATEGORY", "SIMPLE_RACE_CATEGORY1",
  "RAD", "YOST", "histology",
  "baseline_Pneumonitis", "baseline_Hypothyroidism", "baseline_Diarrhea",
  "baseline_AKI", "Charlson_baseline", "adherence_ED_baseline", "Charlson",
  
  # Additional variables used in censoring/weights
  "adherence", "adherence_ED", "adherence_chemo",
  "Pneumonitis", "Hypothyroidism", "Diarrhea", "AKI",
  
  # ID
  "PATIENT_ID",
  
  # Treatment indicator used earlier (if present)
  "CHEMO1"
)

processed_data <- processed_data_raw %>%
  select(any_of(vars_to_keep))

## ---- Replace NA with 0 for Binary/Count Variables ----
vars_replace0 <- c(
  "baseline_Pneumonitis", "baseline_Hypothyroidism", "baseline_Diarrhea",
  "baseline_AKI", "Charlson_baseline", "adherence_ED_baseline",
  "adherence", "adherence_ED", "adherence_chemo",
  "Pneumonitis", "Hypothyroidism", "Diarrhea", "AKI"
)

processed_data <- processed_data %>%
  mutate(across(all_of(intersect(vars_replace0, names(.))),
                ~ if_else(is.na(.x), 0, .x)))

## ---- Define Treatment Category ----
processed_data <- processed_data %>%
  mutate(
    Treatment = case_when(
      adherence_chemo == 1 & adherence == 1 ~ "ICI_Chemo",
      adherence_chemo == 0 & adherence == 1 ~ "ICI_Mono",
      adherence_chemo == 0 & adherence == 0 ~ "No_treatment",
      adherence_chemo == 1 & adherence == 0 ~ "Chemo_only",
      TRUE ~ NA_character_
    )
  )

## ---- Identify CITATION-2 Eligible Patients ----
# Keep patients who received ICI-Chemo in cycles 1–4 and have all four cycles observed.
citation2_ids <- processed_data %>%
  filter(cycle %in% 1:4) %>%
  group_by(PATIENT_ID) %>%
  summarise(
    n_cycle = n_distinct(cycle),
    all_ici_chemo = all(Treatment == "ICI_Chemo"),
    .groups = "drop"
  ) %>%
  filter(n_cycle == 4, all_ici_chemo) %>%
  transmute(PATIENT_ID)

processed_data_citation2 <- processed_data %>%
  filter(PATIENT_ID %in% citation2_ids$PATIENT_ID)

## ---- Re-Index Cycles After Landmark (Start Follow-up After Cycle 4) ----
processed_data_citation2 <- processed_data_citation2 %>%
  mutate(cycle = cycle - 4) %>%
  filter(cycle > 0)

## ---- Exclude Early Death at Follow-up Cycle 1 ----
early_death_ids <- processed_data_citation2 %>%
  filter(cycle == 1, death == 1) %>%
  pull(PATIENT_ID)

processed_data_clean <- processed_data_citation2 %>%
  filter(!PATIENT_ID %in% early_death_ids)

## ---- Create Lagged Outcome (death1) and Time Index ----
# death1 is death in the subsequent cycle interval; remove final row per patient.
processed_long <- processed_data_clean %>%
  group_by(PATIENT_ID) %>%
  arrange(cycle) %>%
  mutate(
    death1 = as.numeric(lead(death, 1)),
    death1 = if_else(is.na(death1), 0, death1),
    time   = row_number() - 1
  ) %>%
  filter(row_number() < n()) %>%
  ungroup()

## ---- Administrative Censoring Flag (Not Used Later, Kept for Completeness) ----
processed_long <- processed_long %>%
  group_by(PATIENT_ID) %>%
  mutate(ad_censor = if_else(row_number() == n() & death1 == 0, 1, 0)) %>%
  ungroup()

## ---- Define Treatment Arm by Cycle 1 After Landmark ----
# chemo_binary: 1 if cycle 1 Treatment == ICI_Chemo; 0 if ICI_Mono.
processed_long <- processed_long %>%
  group_by(PATIENT_ID) %>%
  mutate(
    chemo_binary = case_when(
      any(cycle == 1 & Treatment == "ICI_Chemo", na.rm = TRUE) ~ 1,
      any(cycle == 1 & Treatment == "ICI_Mono",  na.rm = TRUE) ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(chemo_binary))

## ---- Factor / Numeric Casting ----
processed_long <- processed_long %>%
  mutate(
    region                = factor(REGION),
    METRO_CATEGORY        = factor(METRO_CATEGORY),
    SEX_CATEGORY          = factor(SEX_CATEGORY),
    SIMPLE_RACE_CATEGORY1 = factor(SIMPLE_RACE_CATEGORY1),
    histology             = factor(histology),
    adherence_ED          = factor(adherence_ED),
    adherence             = factor(adherence),
    DIAGNOSIS_YEAR1        = factor(DIAGNOSIS_YEAR1),
    adherence_ED_baseline = factor(adherence_ED_baseline),
    
    cycle                 = as.numeric(cycle),
    TTI                   = as.numeric(TTI),
    Charlson_baseline     = as.numeric(Charlson_baseline),
    YOST                  = as.numeric(YOST),
    Charlson              = as.numeric(Charlson)
  )

## ---- Create Grouped Covariates ----
processed_long <- processed_long %>%
  mutate(
    ICI_group = if_else(ICI == "Pembrolizumab", "Pembrolizumab", "Atezolizumab/Nivolumab"),
    ICI_group = factor(ICI_group),
    
    DIAGNOSIS_YEAR1_num = as.numeric(as.character(DIAGNOSIS_YEAR1)),
    DIAGNOSIS_YEAR1_grp = if_else(DIAGNOSIS_YEAR1_num < 2018, "≤2017", "2018–2019"),
    DIAGNOSIS_YEAR1_grp = factor(DIAGNOSIS_YEAR1_grp)
  )

## ---- Truncate Follow-up (Max 48 Cycles = 144 Weeks) ----
processed_long <- processed_long %>% filter(cycle <= 48)

## ---- Split Arms ----
ids_chemo <- processed_long %>%
  filter(cycle == 1, chemo_binary == 1) %>%
  distinct(PATIENT_ID)

df_chemo <- processed_long %>% filter(PATIENT_ID %in% ids_chemo$PATIENT_ID)
df_mono  <- processed_long %>% filter(!PATIENT_ID %in% ids_chemo$PATIENT_ID)

## ---- Per-Protocol with One-Cycle Grace Period ----
add_grace_pp <- function(df, expected_treat) {
  
  dev_info <- df %>%
    group_by(PATIENT_ID) %>%
    summarise(
      exp_treat = expected_treat,
      treat_c2  = Treatment[cycle == 2][1],
      treat_c3  = Treatment[cycle == 3][1],
      
      grace_ok_patient = case_when(
        !is.na(treat_c2) & treat_c2 == exp_treat ~ TRUE,
        !is.na(treat_c2) & treat_c2 == "No_treatment" &
          !is.na(treat_c3) & treat_c3 == exp_treat ~ TRUE,
        TRUE ~ FALSE
      ),
      
      first_dev_cycle = case_when(
        grace_ok_patient ~ Inf,
        !is.na(treat_c2) & treat_c2 != exp_treat & treat_c2 != "No_treatment" ~ 2,
        !is.na(treat_c2) & treat_c2 == "No_treatment" ~ 3,
        TRUE ~ 2
      ),
      .groups = "drop"
    )
  
  df_pp <- df %>%
    left_join(dev_info, by = "PATIENT_ID") %>%
    mutate(
      censor_grace = if_else(
        is.finite(first_dev_cycle) & cycle == first_dev_cycle,
        1L, 0L
      )
    ) %>%
    filter(is.infinite(first_dev_cycle) | cycle <= first_dev_cycle) %>%
    ungroup()
  
  df_pp
}

df_chemo_pp <- add_grace_pp(df_chemo, expected_treat = "ICI_Chemo") %>%
  mutate(arm_name = "4 cycles of ICI-Chemo (one-cycle grace period)")

df_mono_pp <- add_grace_pp(df_mono, expected_treat = "ICI_Mono") %>%
  mutate(arm_name = "ICI-Mono strategy (one-cycle grace period)")

## ---- Main Pipeline Function ----
run_full_pipeline <- function(
    df_arm1,
    df_arm2,
    censor_var,
    out_censor_plot,
    weight_cycles,
    num_formula,
    den_formula,
    out_weight_plot,
    outcome_formula,
    out_boot_rdata,
    out_rd_docx,
    out_surv_plot
) {
  
  ## 1) Censoring flow plot
  flow_df <- bind_rows(df_arm1, df_arm2) %>%
    mutate(status = if_else(.data[[censor_var]] == 1L, 1L, 2L))
  
  plot_counts <- flow_df %>%
    group_by(cycle, status, arm_name) %>%
    summarise(total_count = n(), .groups = "drop")
  
  cens_plot <- ggplot(plot_counts, aes(x = cycle, y = total_count, fill = factor(status))) +
    geom_bar(position = "stack", stat = "identity") +
    xlab("Cycle") +
    ylab("Number of individuals") +
    scale_x_continuous(breaks = seq(0, max(plot_counts$cycle, na.rm = TRUE), by = 4)) +
    scale_fill_manual(
      values = c("1" = "#878787", "2" = "#800000"),
      labels = c("1" = "Censored due to non-adherence", "2" = "Following strategy as planned"),
      name   = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      axis.line = element_line(colour = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    facet_wrap(~ arm_name, scales = "free_y")
  
  ggsave(out_censor_plot, cens_plot, width = 10, height = 5, dpi = 1200, units = "in", bg = "white")
  
  ## 2) Weights for arm 1
  model_data_1 <- df_arm1 %>% filter(cycle %in% weight_cycles)
  
  num_model_1 <- glm(num_formula, family = quasibinomial(), data = model_data_1)
  den_model_1 <- glm(den_formula, family = quasibinomial(), data = model_data_1)
  
  df_arm1_w <- df_arm1 %>%
    mutate(
      num_p = if_else(cycle %in% weight_cycles, predict(num_model_1, newdata = ., type = "response"), 1),
      den_p = if_else(cycle %in% weight_cycles, predict(den_model_1, newdata = ., type = "response"), 1),
      sw_raw = num_p / den_p,
      sw     = if_else(cycle %in% weight_cycles, sw_raw, 1)
    )
  
  lower_1 <- quantile(df_arm1_w$sw_raw, 0.025, na.rm = TRUE)
  upper_1 <- quantile(df_arm1_w$sw_raw, 0.975, na.rm = TRUE)
  
  p1 <- ggplot(df_arm1_w, aes(x = sw)) +
    geom_density(size = 1, color = "black") +
    geom_vline(xintercept = c(lower_1, upper_1), linetype = "dashed", color = "blue", size = 1) +
    coord_cartesian(xlim = c(0.5, 1.5)) +
    labs(title = unique(df_arm1_w$arm_name), x = "Stabilized Weight", y = "Density") +
    theme_bw()
  
  df_arm1_w <- df_arm1_w %>%
    mutate(sw_truncated = pmin(pmax(sw, lower_1), upper_1)) %>%
    group_by(PATIENT_ID) %>%
    arrange(cycle) %>%
    mutate(sw_cum = cumprod(sw_truncated)) %>%
    ungroup()
  
  ## 3) Weights for arm 2
  model_data_2 <- df_arm2 %>% filter(cycle %in% weight_cycles)
  
  num_model_2 <- glm(num_formula, family = quasibinomial(), data = model_data_2)
  den_model_2 <- glm(den_formula, family = quasibinomial(), data = model_data_2)
  
  df_arm2_w <- df_arm2 %>%
    mutate(
      num_p = if_else(cycle %in% weight_cycles, predict(num_model_2, newdata = ., type = "response"), 1),
      den_p = if_else(cycle %in% weight_cycles, predict(den_model_2, newdata = ., type = "response"), 1),
      sw_raw = num_p / den_p,
      sw     = if_else(cycle %in% weight_cycles, sw_raw, 1)
    )
  
  lower_2 <- quantile(df_arm2_w$sw_raw, 0.025, na.rm = TRUE)
  upper_2 <- quantile(df_arm2_w$sw_raw, 0.975, na.rm = TRUE)
  
  p2 <- ggplot(df_arm2_w, aes(x = sw)) +
    geom_density(size = 1, color = "black") +
    geom_vline(xintercept = c(lower_2, upper_2), linetype = "dashed", color = "blue", size = 1) +
    coord_cartesian(xlim = c(0.5, 1.5)) +
    labs(title = unique(df_arm2_w$arm_name), x = "Stabilized Weight", y = "Density") +
    theme_bw()
  
  df_arm2_w <- df_arm2_w %>%
    mutate(sw_truncated = pmin(pmax(sw, lower_2), upper_2)) %>%
    group_by(PATIENT_ID) %>%
    arrange(cycle) %>%
    mutate(sw_cum = cumprod(sw_truncated)) %>%
    ungroup()
  
  ggsave(out_weight_plot, p1 + p2, width = 10, height = 5, dpi = 1200, units = "in", bg = "white")
  
  ## 4) Combine two arms
  df_arm1_w <- df_arm1_w %>% mutate(chemo_binary = 1)
  df_arm2_w <- df_arm2_w %>% mutate(chemo_binary = 0)
  combined_pp <- bind_rows(df_arm1_w, df_arm2_w)
  
  ## 5) Cox models (robust sandwich + cluster)
  cox_144 <- coxph(
    Surv(cycle - 1, cycle, death1) ~
      chemo_binary + ICI + DIAGNOSIS_YEAR1 + age + TTI +
      region + METRO_CATEGORY + SEX_CATEGORY + SIMPLE_RACE_CATEGORY1 + YOST + histology +
      baseline_Pneumonitis + baseline_Hypothyroidism + baseline_Diarrhea +
      baseline_AKI + Charlson_baseline + adherence_ED_baseline,
    data    = combined_pp,
    weights = sw_cum,
    robust  = TRUE,
    cluster = PATIENT_ID
  )
  
  cox_24 <- coxph(
    Surv(cycle - 1, cycle, death1) ~
      chemo_binary + ICI + DIAGNOSIS_YEAR1 + age + TTI +
      region + METRO_CATEGORY + SEX_CATEGORY + SIMPLE_RACE_CATEGORY1 + YOST + histology +
      baseline_Pneumonitis + baseline_Hypothyroidism + baseline_Diarrhea +
      baseline_AKI + Charlson_baseline + adherence_ED_baseline,
    data    = combined_pp %>% filter(cycle <= 8),
    weights = sw_cum,
    robust  = TRUE,
    cluster = PATIENT_ID
  )
  
  cox_48 <- coxph(
    Surv(cycle - 1, cycle, death1) ~
      chemo_binary + DIAGNOSIS_YEAR1 + age + TTI +
      region + METRO_CATEGORY + SEX_CATEGORY + SIMPLE_RACE_CATEGORY1 + YOST + histology +
      baseline_Pneumonitis + baseline_Hypothyroidism + baseline_Diarrhea +
      baseline_AKI + Charlson_baseline + adherence_ED_baseline,
    data    = combined_pp %>% filter(cycle <= 16),
    weights = sw_cum,
    robust  = TRUE,
    cluster = PATIENT_ID
  )
  
  ## 6) Discrete-time hazard model
  discrete_hazard_model <- glm(
    outcome_formula,
    family  = quasibinomial(),
    data    = combined_pp,
    weights = sw_cum
  )
  
  ## 7) Bootstrap standardized survival curves
  baseline <- combined_pp %>%
    group_by(PATIENT_ID) %>%
    filter(cycle == 1) %>%
    ungroup() %>%
    select(
      PATIENT_ID, ICI, DIAGNOSIS_YEAR1, age, TTI,
      region, METRO_CATEGORY, SEX_CATEGORY, SIMPLE_RACE_CATEGORY1,
      YOST, histology,
      baseline_Pneumonitis, baseline_Hypothyroidism, baseline_Diarrhea,
      baseline_AKI, Charlson_baseline, adherence_ED_baseline
    )
  
  categorical_vars <- c(
    "ICI", "DIAGNOSIS_YEAR1", "region", "METRO_CATEGORY", "SEX_CATEGORY",
    "SIMPLE_RACE_CATEGORY1", "histology"
  )
  baseline[categorical_vars] <- lapply(baseline[categorical_vars], factor)
  
  standardization_full <- function(data, indices) {
    base_b <- data[indices, ]
    ids_b  <- base_b$PATIENT_ID
    
    combined_b <- combined_pp %>% filter(PATIENT_ID %in% ids_b)
    
    model_b <- glm(
      outcome_formula,
      family  = quasibinomial(),
      data    = combined_b,
      weights = sw_cum
    )
    
    t_vec <- 1:48
    dfA0  <- data.frame(cycle = t_vec, chemo_binary = 0)
    dfA1  <- data.frame(cycle = t_vec, chemo_binary = 1)
    
    base_b[categorical_vars] <- lapply(base_b[categorical_vars], factor)
    expanded_base <- base_b[rep(seq_len(nrow(base_b)), each = 48), ]
    
    A0 <- cbind(dfA0, expanded_base)
    A1 <- cbind(dfA1, expanded_base)
    
    A0$pred <- predict(model_b, newdata = A0, type = "response")
    A1$pred <- predict(model_b, newdata = A1, type = "response")
    
    n_id <- nrow(base_b)
    A0$id <- rep(seq_len(n_id), each = length(t_vec))
    A1$id <- rep(seq_len(n_id), each = length(t_vec))
    
    A0_surv <- A0 %>% group_by(id) %>% mutate(surv0 = cumprod(1 - pred))
    A1_surv <- A1 %>% group_by(id) %>% mutate(surv1 = cumprod(1 - pred))
    
    s0 <- A0_surv %>% group_by(cycle) %>% summarise(mean_surv0 = mean(surv0), .groups = "drop")
    s1 <- A1_surv %>% group_by(cycle) %>% summarise(mean_surv1 = mean(surv1), .groups = "drop")
    
    merge(s0, s1, by = "cycle")
  }
  
  boot_survival <- function(data, indices) {
    s <- standardization_full(data, indices)
    c(s$mean_surv0, s$mean_surv1)
  }
  
  set.seed(123)
  n_boot <- 300
  boot_results_curves <- boot(data = baseline, statistic = boot_survival, R = n_boot)
  save(boot_results_curves, file = out_boot_rdata)
  
  ## 8) Risk difference table + survival plot
  num_time <- 48
  surv_mat <- boot_results_curves$t
  
  surv0_boot <- surv_mat[, 1:num_time, drop = FALSE]
  surv1_boot <- surv_mat[, (num_time + 1):(2 * num_time), drop = FALSE]
  
  surv0_t0 <- boot_results_curves$t0[1:num_time]
  surv1_t0 <- boot_results_curves$t0[(num_time + 1):(2 * num_time)]
  
  diff_boot <- surv1_boot - surv0_boot
  lower_diff <- apply(diff_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  upper_diff <- apply(diff_boot, 2, quantile, probs = 0.975, na.rm = TRUE)
  
  surv_diff_df <- data.frame(
    cycle      = 1:num_time,
    mean_surv0 = surv0_t0,
    mean_surv1 = surv1_t0,
    mean_diff  = surv1_t0 - surv0_t0,
    lower_diff = lower_diff,
    upper_diff = upper_diff
  )
  
  ft_rd <- flextable(surv_diff_df %>% select(cycle, mean_diff, lower_diff, upper_diff)) %>%
    set_caption("Standardized survival difference over time") %>%
    autofit()
  
  save_as_docx(ft_rd, path = out_rd_docx)
  
  # RD at 144 weeks (cycle 48)
  rd    <- surv_diff_df$mean_diff[num_time]
  rd_lo <- surv_diff_df$lower_diff[num_time]
  rd_hi <- surv_diff_df$upper_diff[num_time]
  
  rd_txt <- sprintf(
    "144-week risk difference: %.1f%% (95%% CI, %.1f%% to %.1f%%)",
    rd * 100, rd_lo * 100, rd_hi * 100
  )
  
  surv_plot_df <- bind_rows(
    data.frame(cycle = 0, mean_surv0 = 1, mean_surv1 = 1),
    surv_diff_df %>% select(cycle, mean_surv0, mean_surv1)
  )
  
  p_surv <- ggplot(surv_plot_df, aes(x = cycle * 3)) +
    geom_line(aes(y = mean_surv0, colour = "Strategy 0"), size = 1) +
    geom_line(aes(y = mean_surv1, colour = "Strategy 1"), size = 1) +
    scale_x_continuous(
      limits = c(0, 144),
      breaks = seq(0, 144, by = 12),
      name   = "Time, weeks"
    ) +
    scale_y_continuous(
      limits = c(0.1, 1),
      breaks = seq(0.1, 1, by = 0.1),
      labels = number_format(accuracy = 0.1),
      name   = "Proportion Alive"
    ) +
    labs(colour = NULL) +
    annotate("text", x = 6, y = 0.15, label = rd_txt, size = 4, hjust = 0) +
    theme_bw(base_family = "Times") +
    theme(
      panel.border       = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "gray85"),
      panel.grid.minor   = element_blank(),
      axis.title         = element_text(size = 11, face = "bold"),
      axis.text          = element_text(size = 10),
      axis.line          = element_line(color = "black"),
      legend.position    = c(0.8, 0.85),
      legend.text        = element_text(size = 10)
    )
  
  ggsave(out_surv_plot, p_surv, width = 6.5, height = 5, dpi = 1200, units = "in", bg = "white")
  
  ## Return results
  invisible(list(
    arm1 = df_arm1_w,
    arm2 = df_arm2_w,
    combined_pp = combined_pp,
    numerator_model_arm1 = summary(num_model_1),
    denominator_model_arm1 = summary(den_model_1),
    numerator_model_arm2 = summary(num_model_2),
    denominator_model_arm2 = summary(den_model_2),
    cox_144week = summary(cox_144),
    cox_24week  = summary(cox_24),
    cox_48week  = summary(cox_48),
    outcome_model_summary = summary(discrete_hazard_model),
    boot_results_curves = boot_results_curves,
    surv_diff_table = surv_diff_df,
    rd_144week_text = rd_txt
  ))
}

## ---- Run Pipeline ----
res_grace <- run_full_pipeline(
  df_arm1         = df_chemo_pp,
  df_arm2         = df_mono_pp,
  censor_var      = "censor_grace",
  out_censor_plot = glue("{save_path}/CITATION-2_cens_plot_grace.png"),
  weight_cycles   = 2:6,
  num_formula     = (censor_grace == 0) ~ cycle + ICI_group + DIAGNOSIS_YEAR1_grp +
    age + TTI + region + METRO_CATEGORY + SEX_CATEGORY + SIMPLE_RACE_CATEGORY1 + YOST + histology,
  den_formula     = (censor_grace == 0) ~ cycle + ICI_group + DIAGNOSIS_YEAR1_grp +
    age + TTI + region + METRO_CATEGORY + SEX_CATEGORY + SIMPLE_RACE_CATEGORY1 + YOST + histology +
    adherence_ED + Pneumonitis + AKI + Charlson,
  out_weight_plot = glue("{save_path}/CITATION-2_weight_grace.png"),
  outcome_formula = death1 ~ ns(cycle, df = 4) +
    chemo_binary + I(ns(cycle, df = 4) * chemo_binary) +
    DIAGNOSIS_YEAR1 + ICI + age + TTI + region +
    METRO_CATEGORY + SEX_CATEGORY + SIMPLE_RACE_CATEGORY1 +
    YOST + histology +
    baseline_Pneumonitis + baseline_Hypothyroidism +
    baseline_Diarrhea + baseline_AKI +
    Charlson_baseline + adherence_ED_baseline,
  out_boot_rdata  = glue("{save_path}/CITATION-2_pp_grace_boot.Rdata"),
  out_rd_docx     = glue("{save_path}/CITATION-2_PP_CI_grace.docx"),
  out_surv_plot   = glue("{save_path}/CITATION-2_PP_grace_survival_curve.png")
)

## ---- Export Model Summaries to DOCX ----
make_model_docx <- function(res_list, path) {
  
  glm_names <- c(
    "numerator_model_arm1",
    "denominator_model_arm1",
    "numerator_model_arm2",
    "denominator_model_arm2"
  )
  
  glm_titles <- c(
    "Numerator model (Arm 1)",
    "Denominator model (Arm 1)",
    "Numerator model (Arm 2)",
    "Denominator model (Arm 2)"
  )
  
  glm_tables <- Map(function(nm, ttl) {
    coefs <- res_list[[nm]]$coefficients
    df <- data.frame(
      Term         = rownames(coefs),
      Estimate     = coefs[, "Estimate"],
      `Std. Error` = coefs[, "Std. Error"],
      row.names    = NULL,
      check.names  = FALSE
    )
    flextable(df) %>% set_caption(ttl) %>% autofit()
  }, glm_names, glm_titles)
  
  cox_names <- c("cox_144week", "cox_24week", "cox_48week")
  cox_titles <- c("Cox model (0–144 weeks)", "Cox model (0–24 weeks)", "Cox model (0–48 weeks)")
  
  cox_tables <- Map(function(nm, ttl) {
    sm <- res_list[[nm]]
    ci <- sm$conf.int
    df <- data.frame(
      Term           = rownames(ci),
      HR             = ci[, "exp(coef)"],
      `Lower 95% CI` = ci[, "lower .95"],
      `Upper 95% CI` = ci[, "upper .95"],
      row.names      = NULL,
      check.names    = FALSE
    )
    flextable(df) %>% set_caption(ttl) %>% autofit()
  }, cox_names, cox_titles)
  
  out_coefs <- res_list$outcome_model_summary$coefficients
  outcome_table <- flextable(
    data.frame(Term = rownames(out_coefs), out_coefs, row.names = NULL, check.names = FALSE)
  ) %>%
    set_caption("Discrete-time outcome model") %>%
    autofit()
  
  args_list <- c(
    setNames(glm_tables, glm_titles),
    setNames(cox_tables, cox_titles),
    setNames(list(outcome_table), "Discrete-time outcome model"),
    list(path = path)
  )
  
  do.call(save_as_docx, args_list)
}

out_model_docx <- glue("{save_path}/CITATION-2_grace_model_summaries.docx")
make_model_docx(res_grace, path = out_model_docx)