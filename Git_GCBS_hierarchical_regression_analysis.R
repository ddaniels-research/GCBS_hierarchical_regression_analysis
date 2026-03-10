# Psychological and Demographic Antecedents of Conspiracy Beliefs:
# A Hierarchical Regression Analysis

# Project:    Generic Conspiracist Beliefs Scale (GCBS) Regression Analysis
# Researcher: Drew Daniels
# Purpose:    Hierarchical regression examining personality, cognitive, and 
#             demographic predictors of conspiracy beliefs, with focus on
#             suppression effects and relative importance analysis

# Theoretical Framework:
#   This analysis approaches conspiracy ideation through a unified, 
#   multi-layered systems framework:
#   1. Interpersonal Layer (Trait Antagonism & Epistemic Mistrust):
#      Individuals low in agreeableness possess a baseline cynicism and 
#      mistrust of others, predisposing them to endorse narratives of 
#      deception and malevolent collusion (H1a).
#   2. Cognitive Layer (Dual-Process Theory & Apophenia):
#      High Openness drives a System 1 heuristic exploration of fringe 
#      ideas. Without the System 2 analytical override provided by higher 
#      cognitive ability (Vocabulary) to reject false patterns (apophenia), 
#      open-mindedness increases conspiracy belief (H1b, H2).
#   3. Environmental Layer (Digital Cohort Effects):
#      Younger cohorts, immersed in decentralized, algorithmically driven 
#      digital ecosystems, exhibit higher baseline exposure and vulnerability 
#      to conspiratorial narratives compared to older cohorts with 
#      crystallized media literacy (H3).

# Data Source: Open Psychometrics GCBS dataset (2016)
# Sample:     Adults who completed online GCBS and agreed to research use

# Key Analytical Features:
#   1. Hierarchical multiple regression (4 blocks)
#   2. Suppression effects detection and interpretation
#   3. Relative importance analysis (dominance analysis)
#   4. Comprehensive assumption checking
#   5. Robustness testing across subgroups


# Ethics (Secondary Data):
#   □ Data collected with informed consent (users opted in for research)
#   □ All participants age 13+ (per codebook)
#   □ No personally identifiable information in dataset
#   □ Data publicly shared for research purposes

#   Research Question: "To what extent do personality and cognitive factors 
# provide incremental predictive value for conspiracy beliefs beyond baseline demographics, 
# and how do cognitive variables suppress or attenuate these relationships?"

#   Hypotheses:
#     H1a (Personality): Lower Agreeableness will be the strongest zero-order
#          personality predictor of conspiracy beliefs.
#     H1b (Suppression): Openness to Experience will demonstrate a positive
#          predictive effect only after controlling for cognitive variables.
#     H2  (Cognitive Mediation): Educational attainment will show a negative
#          association with conspiracy beliefs, significantly attenuated when
#          Vocabulary Knowledge is included.
#     H3  (Age-Linear): Age will demonstrate a significant negative linear
#          relationship with conspiracy beliefs.
#   Analysis Plan:
#     - Hierarchical regression (4 blocks: demographics, personality, cognitive, social)
#     - Alpha = .05, two-tailed tests
#     - Outliers defined as Cook's D > 4/n
#     - Missing data: Listwise deletion (document n at each stage)
#     - Effect sizes: R², ΔR², semi-partial r²
#   Exploratory (clearly labeled):
#     - Subgroup analyses by gender, age groups
#     - Alternative model specifications

# SAMPLE SIZE & POWER (Field et al., 2012, p. 273):
#   With ~15 predictors: minimum N = 50 + 8(15) = 170
#   For medium effect (f² = .15), power = .80: ~150 cases
#   Expected dataset: several thousand (more than adequate)

# ---- STEP 1: PROJECT SETUP & CONFIGURATION ----

# ---- 1.1 Analysis Parameters ----
# Consolidate all decision thresholds here for transparency and easy sensitivity checks

SEED                    <- 42
RAPID_RESPONSE_IQR_MULT <- 1.5     # Tukey lower-fence multiplier for response time
RAPID_RESPONSE_FLOOR    <- 30      # Absolute minimum plausible completion time (seconds)
GCBS_MIN_ITEMS          <- 12      # Require 80% of 15 items for valid score
FAKE_WORD_THRESHOLD     <- 2       # Exclude if ≥ this many fake words endorsed
AGE_RANGE               <- c(13, 100)  # Valid age bounds
SUPPRESSION_MULTIPLIER  <- 1.1     # |β| must exceed |r| × this for suppression
SUPPRESSION_EXPLORE    <- 1.5     # Stricter threshold for exploratory scan
ATTENUATION_THRESHOLD   <- 20      # % reduction considered meaningful for H2
COOKS_D_MULTIPLIER      <- 4       # Cook's D cutoff = multiplier / n
SHAPIRO_SUBSAMPLE_N     <- 5000    # Cap for Shapiro-Wilk (very sensitive with large n)
MULTICOLLINEARITY_R     <- 0.70    # Predictor correlation concern threshold
VIF_MODERATE            <- 5
VIF_SEVERE              <- 10

# ---- 1.2 Helper Functions ----

hr <- function(char = "=", width = 70) {
  cat(strrep(char, width), "\n")
}

log_step <- function(msg) {
  cat(paste0("✓ ", msg, "\n"))
}

format_p <- function(p, digits = 3) {
  format.pval(p, digits = digits)
}

#' Print model fit statistics consistently
print_model_fit <- function(model_summary, prev_summary = NULL) {
  cat("  R²     =", round(model_summary$r.squared, 4), "\n")
  cat("  Adj R² =", round(model_summary$adj.r.squared, 4), "\n")
  if (!is.null(prev_summary)) {
    cat("  ΔR²    =", round(model_summary$r.squared - prev_summary$r.squared, 4), "\n")
  }
  f <- model_summary$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  cat(sprintf("  F(%g, %g) = %.2f, p = %s\n", f[2], f[3], f[1], format_p(p)))
}

#' Print ΔR² F-change test from anova comparison
print_f_change <- function(anova_result) {
  cat(sprintf("  F change(%g, %g) = %.2f, p = %s\n",
              anova_result$Df[2], anova_result$Res.Df[2],
              round(anova_result$F[2], 2),
              format_p(anova_result$`Pr(>F)`[2])))
}

#' Create boxplot of GCBS by a categorical variable
make_boxplot <- function(data, x_var, fill_color, title) {
  ggplot(data, aes(x = .data[[x_var]], y = gcbs_total)) +
    geom_boxplot(fill = fill_color, alpha = 0.7) +
    geom_jitter(alpha = 0.1, width = 0.2) +
    labs(title = title, x = NULL, y = "GCBS Score") +
    theme_minimal()
}

#' Print frequency table for a categorical variable
describe_categorical <- function(data, var, label) {
  cat(paste0("\n", label, ":\n"))
  tbl <- table(data[[var]], useNA = "ifany")
  print(tbl)
  cat("Proportions (%):\n")
  print(round(prop.table(tbl) * 100, 1))
}

#' Build regression formula from DV and predictor vector
build_formula <- function(dv, predictors) {
  reformulate(predictors, response = dv)
}

# Predictor display labels (single source of truth for all tables and plots)
predictor_labels <- c(
  gcbs_total          = "Conspiracy Beliefs (GCBS)",
  age                 = "Age (years)",
  gender              = "Gender (0=M, 1=F)",
  education_num       = "Education",
  urban_num           = "Urban Background",
  english_native_num  = "English Non-Native",
  extraversion        = "Extraversion (TIPI)",
  agreeableness       = "Agreeableness (TIPI)",
  conscientiousness   = "Conscientiousness (TIPI)",
  emotional_stability = "Emotional Stability (TIPI)",
  openness            = "Openness (TIPI)",
  vocab_score         = "Vocabulary Knowledge",
  currently_married   = "Currently Married",
  familysize          = "Family Size"
)

# ---- 1.3 Install and Load Packages ----

install_if_needed <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) install.packages(new_packages, dependencies = TRUE)
}

required_packages <- c(
  # Core data manipulation & workflow
  "here",           # Reproducible file paths
  "tidyverse",      # Data manipulation (dplyr, ggplot2, tidyr, readr, etc.)
  
  # Psychometrics & scale reliability
  "psych",          # Psychological research tools: alpha(), describe()
  
  # Statistical modeling & diagnostics
  "car",            # Companion to Applied Regression: VIF, diagnostic tests
  "lmtest",         # Linear model diagnostic tests (Breusch-Pagan, etc.)
  "lm.beta",        # Standardized regression coefficients
  "sandwich",       # Robust (heteroscedasticity-consistent) standard errors
  "relaimpo",       # Relative importance of predictors in regression models
  
  # Effect sizes
  "effectsize",     # Comprehensive effect size calculations
  
  # Visualization
  "ggcorrplot",     # Correlation matrix plots
  "ggfortify",      # Extends ggplot2 for model diagnostic plots
  "patchwork",      # Combine multiple plots
  
  # Tables & reporting
  "knitr",          # Table formatting
  "kableExtra",     # Enhanced table styling
  "apaTables",      # APA-formatted tables
  
  # Missing data
  "naniar",         # Missing data visualization
  "mice"            # Multiple imputation (available if needed)
)

install_if_needed(required_packages)

library(here)
library(tidyverse)
library(psych)
library(car)
library(lmtest)
library(lm.beta)
library(sandwich)
library(effectsize)
library(relaimpo)
library(ggcorrplot)
library(ggfortify)
library(patchwork)
library(knitr)
library(kableExtra)
library(apaTables)
library(naniar)
library(mice)

# Anchor here() to this script's location
here::i_am("Git_GCBS_hierarchical_regression_analysis.R")

# ---- 1.4 Directory Structure ----

dirs <- c("data/raw", "data/processed", "output/tables", "output/figures",
          "output/models", "scripts", "documentation")
walk(dirs, ~ dir.create(here(.x), recursive = TRUE, showWarnings = FALSE))

log_step("Directory structure created")

# ---- 1.5 Global Options ----

options(scipen = 999, digits = 4)
set.seed(SEED)

# Consistent publication-quality theme
theme_set(theme_minimal(base_size = 12) +
            theme(
              plot.title    = element_text(face = "bold", size = 14),
              plot.subtitle = element_text(size = 11),
              axis.title    = element_text(face = "bold"),
              legend.position = "bottom"
            ))

# ---- 1.6 Session Info (Start) ----

session_start <- sessionInfo()
sink(here("documentation", "session_info_start.txt"))
print(session_start)
sink()

log_step("Session info captured in documentation/session_info_start.txt")


# ---- STEP 2: DATA IMPORT & INITIAL EXPLORATION ----

# ---- 2.1 Import Data ----
# Dataset: GCBS from Open Psychometrics (2016, online administration)

gcbs_raw <- read_csv(
  here("data", "raw", "data.csv"),
  na = c("", "NA", "N/A", "NULL")
)

import_timestamp <- Sys.time()
log_step(paste("Data imported at:", format(import_timestamp)))

# ---- 2.2 Initial Inspection ----

# First 6 rows
print(head(gcbs_raw))

# Dataset Dimensions
cat("Rows:", nrow(gcbs_raw), "| Columns:", ncol(gcbs_raw), "\n")

# Variable Structure
print(str(gcbs_raw, give.attr = FALSE))

# Quick Summary
print(summary(gcbs_raw))

# Check for entirely empty variables
empty_vars <- sapply(gcbs_raw, function(x) all(is.na(x)))
if (any(empty_vars)) {
  cat("\n⚠ WARNING: Completely empty variables:", names(gcbs_raw)[empty_vars], "\n")
} else {
  log_step("No completely empty variables detected")
}

# ---- 2.3 Codebook Reference ----

# GCBS Items (DV):
#   Q1-Q15: 15 items, 1 (definitely not true) to 5 (definitely true)

# Timing (Quality Checks):
#   E1-E15: Time per GCBS item (sec); testelapse: total GCBS time

# TIPI (Big Five, 1–7 scale, two items per dimension, one reverse-scored):
#   E: TIPI1 + (8-TIPI6)     A: (8-TIPI2) + TIPI7     C: TIPI3 + (8-TIPI8)
#   ES: (8-TIPI4) + TIPI9    O: TIPI5 + (8-TIPI10)

# Vocabulary Check:
#   VCL1-VCL16: Word knowledge (1=checked, 0=unchecked)
#   VCL6, VCL9, VCL12 are FAKE WORDS for validity screening

# Demographics:
#   education: 1=<HS, 2=HS, 3=University, 4=Graduate
#   gender: 1=Male, 2=Female, 3=Other
#   age: Years | engnat: 1=Native, 2=Not native
#   urban: 1=Rural, 2=Suburban, 3=Urban
#   voted: 1=Yes, 2=No | married: 1=Never, 2=Currently, 3=Previously
#   familysize: Number of children mother had (including respondent)
#   religion: 1-12 (see codebook) | race: coding issues noted in codebook

# ---- 2.4 Missing Data Overview ----

# Missing Data Overview

missing_summary <- gcbs_raw %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(pct_missing = (n_missing / nrow(gcbs_raw)) * 100) %>%
  arrange(desc(pct_missing))

high_missing <- missing_summary %>% filter(pct_missing > 5)
if (nrow(high_missing) > 0) {
  cat("\n⚠ Variables with >5% missing:\n")
  print(high_missing, n = Inf)
} else {
  log_step("No variables have >5% missing data")
}

# Missing data visualization
missing_plot <- gg_miss_var(gcbs_raw, show_pct = TRUE) +
  labs(title = "Missing Data by Variable",
       subtitle = paste("Dataset:", nrow(gcbs_raw), "cases"),
       x = "Variable", y = "Number Missing")

ggsave(here("output", "figures", "01_missing_data_overview.png"),
       missing_plot, width = 10, height = 8, dpi = 300)
log_step("Missing data visualization saved")

miss_case_result <- miss_case_summary(gcbs_raw)
cat("\nPer-case missingness summary:\n")
print(summary(miss_case_result$n_miss))

# ---- STEP 3: DATA CLEANING & PREPARATION ----

# ---- 3.1 Validity Screening ----

# 3.1.1 Vocabulary check: flag participants endorsing ≥2 fake words (VCL6, VCL9, VCL12)
gcbs_clean <- gcbs_raw %>%
  mutate(
    fake_words_checked = rowSums(pick(VCL6, VCL9, VCL12), na.rm = TRUE),
    invalid_vocabulary = fake_words_checked >= FAKE_WORD_THRESHOLD
  )

n_invalid_vocab <- sum(gcbs_clean$invalid_vocabulary, na.rm = TRUE)

# Vocabulary Validity
cat("Cases endorsing ≥2 fake words:", n_invalid_vocab,
    sprintf("(%.2f%%)\n", n_invalid_vocab / nrow(gcbs_clean) * 100))

# 3.1.2 Age validity: exclude implausible ages or missing
gcbs_clean <- gcbs_clean %>%
  mutate(invalid_age = age < AGE_RANGE[1] | age > AGE_RANGE[2] | is.na(age))

# Age Validity
cat("Invalid age cases:", sum(gcbs_clean$invalid_age, na.rm = TRUE), "\n")

# ---- 3.1.3 Rapid responding (Robust Version) ----
time_median <- median(gcbs_clean$testelapse, na.rm = TRUE)
time_iqr    <- IQR(gcbs_clean$testelapse, na.rm = TRUE)

# Tukey's Outlier Filter (Lower Bound)
time_cutoff  <- time_median - (RAPID_RESPONSE_IQR_MULT * time_iqr)

# Ensure cutoff isn't lower than a "humanly possible" minimum (e.g., 30s)
final_cutoff <- max(time_cutoff, RAPID_RESPONSE_FLOOR)

gcbs_clean <- gcbs_clean %>%
  mutate(rapid_responding = testelapse < final_cutoff | is.na(testelapse))

# Response Time Validity
cat("Median:", round(time_median, 1), "s | IQR:", round(time_iqr, 1), "s\n")
cat("Calculated Bound (Med - 1.5*IQR):", round(time_cutoff, 1), "s\n")
cat("Applied Cutoff (Max of bound or 30s):", round(final_cutoff, 1), "s\n")
cat("Flagged as too fast:", sum(gcbs_clean$rapid_responding, na.rm = TRUE), "\n")

# 3.1.4 Combine all validity flags
gcbs_clean <- gcbs_clean %>%
  mutate(valid_case = !(invalid_vocabulary | invalid_age | rapid_responding))

n_excluded <- sum(!gcbs_clean$valid_case, na.rm = TRUE)

# Overall Validity
cat("Raw N:", nrow(gcbs_raw), "| Excluded:", n_excluded,
    sprintf("(%.2f%%) | Valid:", n_excluded / nrow(gcbs_raw) * 100),
    sum(gcbs_clean$valid_case, na.rm = TRUE), "\n")

# ---- 3.2 Score GCBS (Dependent Variable) ----
# Brotherton et al. (2013): 15 items, mean score, 1–5 scale

gcbs_items <- paste0("Q", 1:15)

# Verify all items exist
missing_items <- gcbs_items[!gcbs_items %in% names(gcbs_clean)]
if (length(missing_items) > 0) {
  stop("Missing GCBS items: ", paste(missing_items, collapse = ", "))
}

# Recode invalid zeros as NA (1–5 scale; zeros are data entry errors)
n_zeros <- sum(gcbs_raw[gcbs_items] == 0, na.rm = TRUE)
gcbs_clean <- gcbs_clean %>%
  mutate(across(all_of(gcbs_items), ~ na_if(.x, 0)))
cat("\nRecoded", n_zeros, "zero values to NA (invalid on 1-5 scale)\n")

# Calculate mean score, requiring ≥80% item completion
gcbs_clean <- gcbs_clean %>%
  mutate(
    gcbs_total  = rowMeans(pick(all_of(gcbs_items)), na.rm = TRUE),
    gcbs_n_items = rowSums(!is.na(pick(all_of(gcbs_items)))),
    gcbs_valid   = gcbs_n_items >= GCBS_MIN_ITEMS
  )

# GCBS Scoring
cat("Valid scores (≥", GCBS_MIN_ITEMS, "items):", sum(gcbs_clean$gcbs_valid), "\n")
cat("Range:", round(range(gcbs_clean$gcbs_total, na.rm = TRUE), 2), "\n")
cat("M =", round(mean(gcbs_clean$gcbs_total, na.rm = TRUE), 2),
    "| SD =", round(sd(gcbs_clean$gcbs_total, na.rm = TRUE), 2), "\n")

# ---- 3.3 Score Big Five (TIPI) ----
# Gosling et al. (2003): 2 items per dimension, one reverse-scored (8 - x)

gcbs_clean <- gcbs_clean %>%
  mutate(
    extraversion        = (TIPI1 + (8 - TIPI6)) / 2,
    agreeableness       = ((8 - TIPI2) + TIPI7) / 2,
    conscientiousness   = (TIPI3 + (8 - TIPI8)) / 2,
    emotional_stability = ((8 - TIPI4) + TIPI9) / 2,
    openness            = (TIPI5 + (8 - TIPI10)) / 2
  )

log_step("TIPI dimensions calculated")

tipi_summary <- gcbs_clean %>%
  summarise(
    across(c(extraversion, agreeableness, conscientiousness, 
             emotional_stability, openness),
           list(mean = ~ mean(., na.rm = TRUE), sd = ~ sd(., na.rm = TRUE)))
  ) %>%
  pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
  separate(stat, into = c("dimension", "statistic"), sep = "_(?=[^_]+$)") %>%
  pivot_wider(names_from = statistic, values_from = value)

# TIPI Summary
print(tipi_summary)

# ---- 3.4 Additional Predictor Variables ----

# 3.4.1 Vocabulary score (real words only, excluding fakes VCL6/9/12)
real_words <- paste0("VCL", c(1:5, 7:8, 10:11, 13:16))

gcbs_clean <- gcbs_clean %>%
  mutate(
    vocab_score       = rowSums(pick(all_of(real_words)), na.rm = TRUE),
    vocab_total_items = rowSums(!is.na(pick(all_of(real_words))))
  )

# Vocabulary Score
cat("Real words:", length(real_words), "| Range:",
    range(gcbs_clean$vocab_score, na.rm = TRUE), "\n")
cat("M =", round(mean(gcbs_clean$vocab_score, na.rm = TRUE), 2), "\n")

# 3.4.2 Recode categorical variables
gcbs_clean <- gcbs_clean %>%
  mutate(
    education_factor = factor(education, levels = 1:4, ordered = TRUE,
                              labels = c("Less than HS", "High School", 
                                         "University", "Graduate")),
    education_num = as.numeric(education),
    
    gender_factor = factor(gender, levels = 1:3,
                           labels = c("Male", "Female", "Other")),
    
    urban_factor = factor(urban, levels = 1:3, ordered = TRUE,
                          labels = c("Rural", "Suburban", "Urban")),
    urban_num = as.numeric(urban),
    
    english_native = factor(engnat, levels = 1:2, labels = c("Yes", "No")),
    
    voted_factor = factor(voted, levels = 1:2, labels = c("Yes", "No")),
    
    married_factor = factor(married, levels = 1:3,
                            labels = c("Never", "Currently", "Previously")),
    currently_married = ifelse(married == 2, 1, 0)
  )

log_step("Categorical variables recoded")

# 3.4.3 Collapse religion (12 categories → 4)
gcbs_clean <- gcbs_clean %>%
  mutate(
    religion_detailed = factor(religion, levels = 1:12,
                               labels = c("Agnostic", "Atheist", "Buddhist",
                                          "Catholic", "Mormon", "Protestant",
                                          "Christian Other", "Hindu", "Jewish",
                                          "Muslim", "Sikh", "Other")),
    religion_collapsed = case_when(
      religion %in% c(1, 2)          ~ "Non-religious",
      religion %in% c(4, 5, 6, 7)    ~ "Christian",
      religion %in% c(3, 8, 9, 10, 11) ~ "Other religion",
      religion == 12                  ~ "Other",
      TRUE                            ~ NA_character_
    ) %>% factor(levels = c("Non-religious", "Christian", "Other religion", "Other"))
  )

cat("\nReligion (collapsed):\n")
print(table(gcbs_clean$religion_collapsed, useNA = "ifany"))

# ---- 3.5 Create Analysis Dataset ----

gcbs_analysis <- gcbs_clean %>%
  filter(valid_case == TRUE, gcbs_valid == TRUE) %>%
  dplyr::select(
    gcbs_total,
    extraversion, agreeableness, conscientiousness, emotional_stability, openness,
    vocab_score,
    age, gender, gender_factor,
    education_num, education_factor,
    urban_num, urban_factor,
    english_native,
    currently_married, married_factor,
    familysize,
    voted_factor,
    religion_collapsed,
    fake_words_checked, testelapse
  )

n_final <- nrow(gcbs_analysis)

# Final Analysis Sample
cat("N =", n_final, sprintf("(%.1f%% of raw)\n", n_final / nrow(gcbs_raw) * 100))

write_csv(gcbs_analysis, here("data", "processed", "gcbs_analysis_ready.csv"))
log_step("Analysis-ready dataset saved")

# ---- 3.6 Document Cleaning Decisions ----

cleaning_summary <- tibble(
  Step = c("Raw data imported",
           "Removed: Failed vocabulary check",
           "Removed: Invalid age",
           "Removed: Rapid responding",
           "Removed: Insufficient GCBS items",
           "Final analysis sample"),
  N_Affected = c(NA,
                 sum(gcbs_clean$invalid_vocabulary, na.rm = TRUE),
                 sum(gcbs_clean$invalid_age, na.rm = TRUE),
                 sum(gcbs_clean$rapid_responding, na.rm = TRUE),
                 sum(!gcbs_clean$gcbs_valid[gcbs_clean$valid_case], na.rm = TRUE),
                 NA),
  N_Remaining = c(nrow(gcbs_raw),
                  nrow(gcbs_raw) - sum(gcbs_clean$invalid_vocabulary, na.rm = TRUE),
                  nrow(gcbs_raw) - sum(gcbs_clean$invalid_age, na.rm = TRUE),
                  nrow(gcbs_raw) - sum(gcbs_clean$rapid_responding, na.rm = TRUE),
                  sum(gcbs_clean$gcbs_valid[gcbs_clean$valid_case], na.rm = TRUE),
                  n_final),
  Note = c("", "≥2 fake words endorsed", "Outside 13–100 or missing",
           "< Tukey lower fence (Med−1.5×IQR) or 30s floor", "< 12 of 15 items",
           "Exclusions applied simultaneously; rows are not sequential")
)

write_csv(cleaning_summary, here("output", "tables", "01_data_cleaning_flowchart.csv"))
log_step("Cleaning summary saved")

# ---- STEP 4: DESCRIPTIVE STATISTICS & VISUALIZATION ----

# DESCRIPTIVE STATISTICS

# ---- 4.1 DV: GCBS ----

# DEPENDENT VARIABLE: Conspiracy Beliefs (GCBS)

gcbs_desc <- gcbs_analysis %>%
  summarise(
    n        = sum(!is.na(gcbs_total)),
    Mean     = mean(gcbs_total, na.rm = TRUE),
    SD       = sd(gcbs_total, na.rm = TRUE),
    Median   = median(gcbs_total, na.rm = TRUE),
    Min      = min(gcbs_total, na.rm = TRUE),
    Max      = max(gcbs_total, na.rm = TRUE),
    Skewness = psych::skew(gcbs_total, na.rm = TRUE),
    Kurtosis = psych::kurtosi(gcbs_total, na.rm = TRUE),
    SE       = SD / sqrt(n),
    CI_lower = Mean - 1.96 * SE,
    CI_upper = Mean + 1.96 * SE
  )
print(gcbs_desc)

cat("Skewness:", ifelse(abs(gcbs_desc$Skewness) < 0.5, "Minimal",
                        ifelse(abs(gcbs_desc$Skewness) < 1, "Moderate", "Severe")), "\n")
cat("Kurtosis:", ifelse(abs(gcbs_desc$Kurtosis) < 0.5, "Minimal",
                        ifelse(abs(gcbs_desc$Kurtosis) < 1, "Moderate", "Severe")), "\n")

# ---- 4.2 Continuous Predictors ----

continuous_vars <- c("extraversion", "agreeableness", "conscientiousness",
                     "emotional_stability", "openness", "vocab_score", "age")

# CONTINUOUS PREDICTORS

continuous_desc <- gcbs_analysis %>%
  dplyr::select(all_of(continuous_vars)) %>%
  psych::describe() %>%
  as.data.frame() %>%
  dplyr::select(n, mean, sd, median, min, max, skew, kurtosis, se) %>%
  round(2)
print(continuous_desc)

# ---- 4.3 Categorical Predictors ----

# CATEGORICAL PREDICTORS

describe_categorical(gcbs_analysis, "gender_factor",     "Gender")
describe_categorical(gcbs_analysis, "education_factor",  "Education")
describe_categorical(gcbs_analysis, "urban_factor",      "Urban/Rural Background")
describe_categorical(gcbs_analysis, "english_native",    "English Native Speaker")
describe_categorical(gcbs_analysis, "voted_factor",      "Voted in Past Year")
describe_categorical(gcbs_analysis, "religion_collapsed", "Religion (Collapsed)")

# ---- 4.4 Distribution Plots ----

# DV distribution with normal overlay
gcbs_dist_plot <- ggplot(gcbs_analysis, aes(x = gcbs_total)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 fill = "steelblue", alpha = 0.7, color = "white") +
  geom_density(color = "darkblue", linewidth = 1) +
  stat_function(fun = dnorm,
                args = list(mean = mean(gcbs_analysis$gcbs_total, na.rm = TRUE),
                            sd = sd(gcbs_analysis$gcbs_total, na.rm = TRUE)),
                color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean(gcbs_total, na.rm = TRUE)),
             color = "darkgreen", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = median(gcbs_total, na.rm = TRUE)),
             color = "orange", linetype = "dotted", linewidth = 1) +
  labs(title = "Distribution of Conspiracy Beliefs (GCBS)",
       subtitle = sprintf("N = %d | M = %.2f | SD = %.2f", n_final,
                          mean(gcbs_analysis$gcbs_total, na.rm = TRUE),
                          sd(gcbs_analysis$gcbs_total, na.rm = TRUE)),
       x = "GCBS Total Score (1–5)", y = "Density",
       caption = "Green dashed = mean | Orange dotted = median | Red dashed = normal") +
  theme_minimal()

gcbs_qq_plot <- ggplot(gcbs_analysis, aes(sample = gcbs_total)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot: GCBS Normality", x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_minimal()

gcbs_dist_combined <- gcbs_dist_plot / gcbs_qq_plot +
  plot_annotation(title = "Conspiracy Beliefs: Distribution Analysis",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave(here("output", "figures", "02_gcbs_distribution.png"),
       gcbs_dist_combined, width = 10, height = 10, dpi = 300)
log_step("GCBS distribution plots saved")

# Personality distributions
personality_vars <- c("extraversion", "agreeableness", "conscientiousness",
                      "emotional_stability", "openness")

personality_plots <- map(personality_vars, function(var) {
  ggplot(gcbs_analysis, aes(x = .data[[var]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "skyblue", 
                   color = "white", alpha = 0.7) +
    geom_density(color = "darkblue", linewidth = 1) +
    labs(title = str_to_title(str_replace_all(var, "_", " ")),
         x = "Score (1–7)", y = "Density") +
    theme_minimal(base_size = 10)
})

personality_combined <- wrap_plots(personality_plots, ncol = 2) +
  plot_annotation(title = "Big Five Personality Distributions (TIPI)",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(here("output", "figures", "03_personality_distributions.png"),
       personality_combined, width = 12, height = 10, dpi = 300)
log_step("Personality distribution plots saved")

# ---- 4.5 Bivariate Relationships ----

# Scatterplots: continuous predictors vs GCBS
scatter_plots <- map(continuous_vars, function(var) {
  ggplot(gcbs_analysis, aes(x = .data[[var]], y = gcbs_total)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", color = "blue", se = TRUE) +
    geom_smooth(method = "loess", color = "red", se = FALSE, linetype = "dashed") +
    labs(title = str_to_title(str_replace_all(var, "_", " ")),
         x = str_to_title(str_replace_all(var, "_", " ")), y = "GCBS Score",
         caption = "Blue = linear | Red = loess") +
    theme_minimal(base_size = 9)
})

scatter_combined <- wrap_plots(scatter_plots, ncol = 3) +
  plot_annotation(title = "Bivariate Relationships: Predictors vs Conspiracy Beliefs",
                  subtitle = "Checking linearity assumptions",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(here("output", "figures", "04_bivariate_scatterplots.png"),
       scatter_combined, width = 14, height = 10, dpi = 300)
log_step("Bivariate relationship plots saved")

# Boxplots: categorical predictors vs GCBS
edu_box     <- make_boxplot(gcbs_analysis, "education_factor",  "lightblue",  "GCBS by Education")
gender_box  <- make_boxplot(gcbs_analysis, "gender_factor",     "lightgreen", "GCBS by Gender")
urban_box   <- make_boxplot(gcbs_analysis, "urban_factor",      "lightyellow","GCBS by Urban/Rural")
religion_box <- make_boxplot(gcbs_analysis, "religion_collapsed","lightpink",  "GCBS by Religion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

categorical_combined <- (edu_box + gender_box) / (urban_box + religion_box) +
  plot_annotation(title = "GCBS by Categorical Predictors",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(here("output", "figures", "05_categorical_predictors_boxplots.png"),
       categorical_combined, width = 12, height = 10, dpi = 300)
log_step("Categorical predictor plots saved")

# ---- 4.6 Publication-Quality Descriptive Table ----

desc_table <- gcbs_analysis %>%
  dplyr::select(gcbs_total, all_of(continuous_vars)) %>%
  psych::describe() %>%
  as.data.frame() %>%
  dplyr::select(n, mean, sd, min, max) %>%
  rownames_to_column("Variable") %>%
  mutate(
    Variable = predictor_labels[Variable],
    across(where(is.numeric), ~ round(., 2))
  )

write_csv(desc_table, here("output", "tables", "02_descriptive_statistics.csv"))

kable(desc_table, caption = "Table 1. Descriptive Statistics for Key Variables") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  save_kable(here("output", "tables", "02_descriptive_statistics.html"))

log_step("Descriptive statistics table saved")

# ---- STEP 5: SCALE RELIABILITY & PSYCHOMETRIC VALIDATION ----

# GCBS RELIABILITY ANALYSIS

# ---- 5.1 Internal Consistency ----

gcbs_items_data <- gcbs_clean %>%
  filter(valid_case == TRUE, gcbs_valid == TRUE) %>%
  dplyr::select(all_of(gcbs_items))

# Cronbach's alpha
alpha_result <- psych::alpha(gcbs_items_data, check.keys = TRUE)

cat("\n--- Cronbach's Alpha ---\n")
cat("Raw α:", round(alpha_result$total$raw_alpha, 3), "\n")
cat("Standardized α:", round(alpha_result$total$std.alpha, 3), "\n")

# McDonald's omega (often more appropriate than alpha)
omega_result <- psych::omega(gcbs_items_data, nfactors = 1, plot = FALSE)

cat("\n--- McDonald's Omega ---\n")
cat("ω total:", round(omega_result$omega.tot, 3), "\n")
cat("ω hierarchical:", round(omega_result$omega_h, 3), "\n")

alpha_val <- alpha_result$total$std.alpha
cat("\nReliability is",
    ifelse(alpha_val >= 0.90, "Excellent (α ≥ .90)",
           ifelse(alpha_val >= 0.80, "Good (α ≥ .80)",
                  ifelse(alpha_val >= 0.70, "Acceptable (α ≥ .70)",
                         "Questionable (α < .70)"))), "\n")

# ---- 5.2 Item-Total Correlations ----

# Item-Total Correlations
# Items with r < .30 may be problematic

item_total <- alpha_result$item.stats %>%
  dplyr::select(raw.r, std.r, `r.cor`, `r.drop`) %>%
  round(3)
print(item_total)

problematic <- which(item_total$r.cor < 0.30)
if (length(problematic) > 0) {
  cat("\n Items with corrected r < .30:", names(problematic), "\n")
} else {
  log_step("All items show adequate item-total correlations (r ≥ .30)")
}

# ---- 5.3 Alpha if Item Deleted ----

# Alpha if Item Deleted

alpha_if_deleted <- alpha_result$alpha.drop %>%
  dplyr::select(`raw_alpha`, `std.alpha`) %>%
  round(3)
print(alpha_if_deleted)

max_improve <- max(alpha_if_deleted$std.alpha) - alpha_result$total$std.alpha
if (max_improve > 0.02) {
  cat("\n Consider removing:",
      rownames(alpha_if_deleted)[which.max(alpha_if_deleted$std.alpha)],
      "(Δα =", round(max_improve, 3), ")\n")
} else {
  log_step("No single item removal would meaningfully improve reliability")
}

# ---- 5.4 Inter-Item Correlations ----

item_cors <- cor(gcbs_items_data, use = "pairwise.complete.obs")
mean_inter_item <- mean(item_cors[lower.tri(item_cors)])

cat("\n--- Inter-Item Correlations ---\n")
cat("Mean inter-item r:", round(mean_inter_item, 3),
    "-", ifelse(mean_inter_item < 0.15, "Too low",
                ifelse(mean_inter_item > 0.50, "Too high (possible redundancy)",
                       "Appropriate range (.15–.50)")), "\n")

# ---- 5.5 Dimensionality Check ----

# Dimensionality (Parallel Analysis)

parallel <- psych::fa.parallel(gcbs_items_data, fm = "ml", fa = "both",
                               n.iter = 100, plot = FALSE)

cat("Suggested factors:", parallel$nfact, "| Components:", parallel$ncomp, "\n")
if (parallel$nfact == 1) {
  log_step("GCBS appears unidimensional")
} else {
  cat("Multiple factors suggested — consider examining subscales\n")
}

# ---- 5.6 Save Reliability ----

reliability_summary <- tibble(
  Measure = "GCBS (15 items)",
  N = nrow(gcbs_items_data),
  `Cronbach's α` = round(alpha_result$total$std.alpha, 3),
  `McDonald's ω` = round(omega_result$omega.tot, 3),
  `Mean Inter-Item r` = round(mean_inter_item, 3),
  Range = paste0(round(min(gcbs_items_data, na.rm = TRUE), 1), "–",
                 round(max(gcbs_items_data, na.rm = TRUE), 1))
)

write_csv(reliability_summary, here("output", "tables", "03_gcbs_reliability.csv"))
log_step("Reliability summary saved")

# ---- STEP 6: CORRELATION ANALYSIS & MULTICOLLINEARITY CHECK ----

# Goals before regression:
# 1. Examine zero-order correlations with GCBS
# 2. Check multicollinearity among predictors
# 3. Identify potential suppression effects (|r with DV| ≈ 0 but correlated with predictors)

# CORRELATION ANALYSIS

# ---- 6.1 Compute Correlation Matrix ----

cor_vars <- c("gcbs_total", "extraversion", "agreeableness", "conscientiousness",
              "emotional_stability", "openness", "vocab_score", "age",
              "education_num", "urban_num", "familysize")

# Note: gender, english_native_num, and currently_married are binary in
# regression_data but haven't been recoded yet in gcbs_analysis at this point.
# Point-biserial correlations for these are computed during regression (Step 7).

cor_data <- gcbs_analysis %>%
  dplyr::select(all_of(cor_vars)) %>%
  drop_na()

cat("N for correlations:", nrow(cor_data), "\n")

cor_matrix       <- cor(cor_data, use = "complete.obs", method = "pearson")
cor_test_results <- psych::corr.test(cor_data, use = "complete", method = "pearson")
p_matrix         <- cor_test_results$p

# ---- 6.2 Zero-Order Correlations with GCBS ----

# ZERO-ORDER CORRELATIONS WITH GCBS

gcbs_cors <- cor_matrix["gcbs_total", ] %>%
  enframe(name = "Variable", value = "r") %>%
  filter(Variable != "gcbs_total") %>%
  mutate(
    p = p_matrix["gcbs_total", Variable],
    sig = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ ""),
    r_display = paste0(round(r, 3), sig)
  ) %>%
  arrange(desc(abs(r)))

print(gcbs_cors)

# Effect sizes (Cohen, 1988): Small |r| .10–.29 | Medium .30–.49 | Large ≥ .50

# ---- 6.3 Multicollinearity Among Predictors ----

# MULTICOLLINEARITY ASSESSMENT

predictor_vars <- setdiff(cor_vars, "gcbs_total")
predictor_cors <- cor_matrix[predictor_vars, predictor_vars]

high_cors <- which(abs(predictor_cors) > MULTICOLLINEARITY_R & 
                     predictor_cors != 1, arr.ind = TRUE)

if (nrow(high_cors) > 0) {
  cat("HIGH CORRELATIONS (|r| >", MULTICOLLINEARITY_R, "):\n")
  high_cor_pairs <- data.frame(
    Var1 = rownames(predictor_cors)[high_cors[, 1]],
    Var2 = colnames(predictor_cors)[high_cors[, 2]],
    r = predictor_cors[high_cors]
  ) %>%
    distinct() %>%
    filter(as.numeric(factor(Var1)) < as.numeric(factor(Var2)))
  print(high_cor_pairs)
} else {
  log_step(paste("No problematic multicollinearity (all |r| <", MULTICOLLINEARITY_R, ")"))
}

# ---- 6.4 H1a Test: Agreeableness as Strongest Zero-Order Personality Predictor ----

# H1a TEST: Agreeableness as Strongest Zero-Order Personality Predictor

personality_cors <- gcbs_cors %>%
  filter(Variable %in% personality_vars) %>%
  mutate(abs_r = abs(r)) %>%
  arrange(desc(abs_r))

# Personality |r| ranking with GCBS:
for (i in 1:nrow(personality_cors)) {
  cat(sprintf("%d. %s: r = %.4f %s\n",
              i, personality_cors$Variable[i],
              personality_cors$r[i], personality_cors$sig[i]))
}

strongest_personality <- personality_cors$Variable[1]
agreeableness_r    <- personality_cors$r[personality_cors$Variable == "agreeableness"]
agreeableness_rank <- which(personality_cors$Variable == "agreeableness")
h1a_supported      <- strongest_personality == "agreeableness"

cat("\n--- H1a RESULT ---\n")
cat("Strongest:", strongest_personality, "| Agreeableness rank:", agreeableness_rank, "\n")
cat("H1a Supported:", ifelse(h1a_supported, "YES ✓",
                             paste0("NO — ", strongest_personality, " is stronger")), "\n")

write_csv(personality_cors, here("output", "tables", "06a_personality_zero_order_correlations.csv"))

# ---- 6.5 Correlation Matrix Visualization ----

cor_plot <- ggcorrplot(
  cor_matrix, hc.order = TRUE, type = "lower", lab = TRUE, lab_size = 3,
  p.mat = p_matrix, sig.level = 0.05, insig = "blank",
  colors = c("#6D9EC1", "white", "#E46726"),
  title = "Correlation Matrix: GCBS and Predictors",
  ggtheme = theme_minimal()
) +
  labs(subtitle = "Only significant correlations shown (p < .05)")

ggsave(here("output", "figures", "06_correlation_matrix.png"),
       cor_plot, width = 12, height = 10, dpi = 300)
log_step("Correlation matrix plot saved")

# ---- 6.6 Identify Potential Suppressors ----

# POTENTIAL SUPPRESSION EFFECTS
# Variables with |r| < .15 with GCBS but correlated with predictors

suppressor_candidates <- gcbs_cors %>% filter(abs(r) < 0.15)

if (nrow(suppressor_candidates) > 0) {
  cat("Candidates:\n")
  for (var in suppressor_candidates$Variable) {
    other_cors <- cor_matrix[var, predictor_vars[predictor_vars != var]]
    max_cor <- max(abs(other_cors))
    cat("  ", var, ": max |r| with other predictors =", round(max_cor, 3))
    if (max_cor > 0.30) cat(" → Potential suppressor!")
    cat("\n")
  }
} else {
  log_step("No obvious suppressor candidates (all |r| > .15 with GCBS)")
}

# ---- 6.7 Save Correlation Tables ----

write_csv(
  as.data.frame(cor_matrix) %>% rownames_to_column("Variable") %>%
    mutate(across(where(is.numeric), ~ round(., 3))),
  here("output", "tables", "04_correlation_matrix_full.csv")
)

write_csv(gcbs_cors, here("output", "tables", "05_gcbs_zero_order_correlations.csv"))
log_step("Correlation tables saved")

# ---- STEP 7: HIERARCHICAL MULTIPLE REGRESSION — MAIN ANALYSIS ----

# Blocks reflect theoretical groupings:
# 1. Demographics: baseline control variables
# 2. + Personality: Big Five dimensions
# 3. + Cognitive: vocabulary knowledge (proxy for cognitive ability)
# 4. + Social: marital status, family size

# HIERARCHICAL MULTIPLE REGRESSION

# 7.0.1 Gender category examination for analytic decisions
# Gender is nominal (1=Male, 2=Female, 3=Other); treating it as continuous is invalid.
# Examine category sizes to determine coding strategy.

# Gender Distribution for Regression Coding

gender_breakdown <- gcbs_analysis %>%
  count(gender, gender_factor) %>%
  mutate(pct = round(n / sum(n) * 100, 2))

print(gender_breakdown)

n_other_gender <- gender_breakdown %>% filter(gender == 3) %>% pull(n)
pct_other_gender <- gender_breakdown %>% filter(gender == 3) %>% pull(pct)

cat(sprintf("\nGender = Other: n = %d (%.2f%% of analysis sample)\n",
            n_other_gender, pct_other_gender))
cat("Decision: Exclude Gender = Other from regression due to insufficient cell size\n")
cat("for stable coefficient estimation. Retained as 0 = Male, 1 = Female.\n")
cat("This exclusion is documented in the data cleaning flowchart.\n\n")

# ---- 7.1 Prepare Complete-Cases Regression Dataset ----
# Enforce identical samples across all models for valid ΔR² comparisons

# Define predictor blocks
block_demographics <- c("age", "gender", "education_num", "urban_num", "english_native_num")
block_personality  <- c("extraversion", "agreeableness", "conscientiousness",
                        "emotional_stability", "openness")
block_cognitive    <- c("vocab_score")
block_social       <- c("currently_married", "familysize")

# Cumulative predictor sets for each block
preds_m1 <- block_demographics
preds_m2 <- c(preds_m1, block_personality)
preds_m3 <- c(preds_m2, block_cognitive)
preds_m4 <- c(preds_m3, block_social)

regression_data <- gcbs_analysis %>%
  filter(gender %in% c(1, 2)) %>%
  dplyr::select(gcbs_total, age, gender, education_num, urban_num, english_native,
                all_of(c(block_personality, block_cognitive, block_social))) %>%
  mutate(
    english_native_num = as.numeric(english_native) - 1,
    gender = gender - 1
  ) %>%
  dplyr::select(-english_native) %>%
  drop_na()

n_regression <- nrow(regression_data)

# Recompute zero-order correlations on the exact regression sample
# so that r values used in suppression/hypothesis tests match the β sample
regression_cor_vars <- c("gcbs_total", preds_m4)
regression_cor_matrix <- cor(regression_data[regression_cor_vars], use = "complete.obs")
min_n <- 50 + 8 * length(c(block_demographics, block_personality, block_cognitive, block_social))

cat("Regression N:", n_regression, "\n")
if (n_regression >= min_n) {
  log_step(paste("Sample size adequate (n =", n_regression, "≥", min_n, ")"))
} else {
  cat("Sample size may be inadequate (n =", n_regression, "<", min_n, ")\n")
}

# ---- 7.2 Build Hierarchical Models ----

# Building Hierarchical Models

# Model 1: Demographics
cat("Model 1: Demographics\n")
model1 <- lm(build_formula("gcbs_total", preds_m1), data = regression_data)
model1_summary <- summary(model1)
model1_std <- lm.beta(model1)
print_model_fit(model1_summary)

# Model 2: + Personality
cat("Model 2: + Personality (Big Five)\n")
model2 <- lm(build_formula("gcbs_total", preds_m2), data = regression_data)
model2_summary <- summary(model2)
model2_std <- lm.beta(model2)
print_model_fit(model2_summary, model1_summary)
anova_1_2 <- anova(model1, model2)
print_f_change(anova_1_2)

# Model 3: + Cognitive
cat("Model 3: + Cognitive (Vocabulary)\n")
model3 <- lm(build_formula("gcbs_total", preds_m3), data = regression_data)
model3_summary <- summary(model3)
model3_std <- lm.beta(model3)
print_model_fit(model3_summary, model2_summary)
anova_2_3 <- anova(model2, model3)
print_f_change(anova_2_3)

# ---- 7.2.1 H1b: Openness Suppression Effect ----
# Classical suppression: |β| > |r| after controlling for cognitive variables

# H1b TEST: Openness Suppression Effect (Classical Suppression)

openness_r <- regression_cor_matrix["gcbs_total", "openness"]
openness_beta_model2 <- coef(model2_std)["openness"]
openness_p_model2    <- model2_summary$coefficients["openness", "Pr(>|t|)"]
openness_beta_model3 <- coef(model3_std)["openness"]
openness_p_model3    <- model3_summary$coefficients["openness", "Pr(>|t|)"]

# Openness Coefficient Tracking:
cat(sprintf("  Zero-order r:          %.4f\n", openness_r))
cat(sprintf("  Model 2 β (pre-cog):   %.4f (p = %s)\n",
            openness_beta_model2, format_p(openness_p_model2)))
cat(sprintf("  Model 3 β (post-cog):  %.4f (p = %s)\n\n",
            openness_beta_model3, format_p(openness_p_model3)))

beta_increase <- openness_beta_model3 - openness_beta_model2
beta_vs_r     <- abs(openness_beta_model3) - abs(openness_r)

# Suppression Indicators:
cat(sprintf("  β change (M2 → M3):   %.4f\n", beta_increase))
cat(sprintf("  |β| - |r|:            %.4f\n\n", beta_vs_r))

suppression_detected  <- abs(openness_beta_model3) > abs(openness_r) * SUPPRESSION_MULTIPLIER
positive_after_control <- openness_beta_model3 > 0 & openness_p_model3 < 0.05
enhancement_detected  <- abs(openness_beta_model3) > abs(openness_beta_model2) * SUPPRESSION_MULTIPLIER

# H1b RESULT ---\n")
cat("Classical suppression (|β| > |r|):", suppression_detected, "\n")
cat("Enhancement when vocab added:", enhancement_detected, "\n")
cat("Positive & significant after control:", positive_after_control, "\n")
cat("H1b Supported:", ifelse(suppression_detected & positive_after_control,
                             "YES ✓ — Suppression effect detected",
                             "NO — Suppression criteria not met"), "\n")

write_csv(
  tibble(
    Statistic = c("Zero-Order r", "Model 2 β", "Model 2 p",
                  "Model 3 β", "Model 3 p", "β change", "|β| - |r|"),
    Value = c(round(openness_r, 4), round(openness_beta_model2, 4),
              format_p(openness_p_model2, 4),
              round(openness_beta_model3, 4), format_p(openness_p_model3, 4),
              round(beta_increase, 4), round(beta_vs_r, 4))
  ),
  here("output", "tables", "07a_openness_suppression_analysis.csv")
)
log_step("Openness suppression analysis saved")

# ---- 7.2.2 H2: Education Attenuation by Cognitive Ability ----

# H2 TEST: Education Effect Attenuation by Cognitive Ability

edu_beta_model1 <- coef(model1_std)["education_num"]
edu_beta_model2 <- coef(model2_std)["education_num"]
edu_beta_model3 <- coef(model3_std)["education_num"]

edu_p_model1 <- model1_summary$coefficients["education_num", "Pr(>|t|)"]
edu_p_model2 <- model2_summary$coefficients["education_num", "Pr(>|t|)"]
edu_p_model3 <- model3_summary$coefficients["education_num", "Pr(>|t|)"]

edu_b_model1 <- coef(model1)["education_num"]
edu_b_model2 <- coef(model2)["education_num"]
edu_b_model3 <- coef(model3)["education_num"]

# Education β Across Models:
cat(sprintf("  Model 1 (demographics): β = %.4f, b = %.4f, p = %s\n",
            edu_beta_model1, edu_b_model1, format_p(edu_p_model1)))
cat(sprintf("  Model 2 (+ personality): β = %.4f, b = %.4f, p = %s\n",
            edu_beta_model2, edu_b_model2, format_p(edu_p_model2)))
cat(sprintf("  Model 3 (+ cognitive):   β = %.4f, b = %.4f, p = %s\n\n",
            edu_beta_model3, edu_b_model3, format_p(edu_p_model3)))

attenuation_pct <- ((abs(edu_beta_model2) - abs(edu_beta_model3)) /
                      abs(edu_beta_model2)) * 100
attenuation_from_baseline <- ((abs(edu_beta_model1) - abs(edu_beta_model3)) /
                                abs(edu_beta_model1)) * 100

# Attenuation:
cat(sprintf("  Model 2 → Model 3: %.1f%%\n", attenuation_pct))
cat(sprintf("  Model 1 → Model 3: %.1f%%\n\n", attenuation_from_baseline))

h2_negative_effect <- edu_beta_model1 < 0 & edu_p_model1 < 0.05
h2_attenuated      <- attenuation_pct > ATTENUATION_THRESHOLD
h2_vocab_significant <- model3_summary$coefficients["vocab_score", "Pr(>|t|)"] < 0.05

# H2 RESULT
cat("Negative education effect in Model 1:", h2_negative_effect, "\n")
cat("Vocabulary significantly predicts GCBS:", h2_vocab_significant, "\n")
cat(sprintf("Meaningful attenuation (>%d%%): %s (%.1f%%)\n",
            ATTENUATION_THRESHOLD, h2_attenuated, attenuation_pct))

if (h2_negative_effect & h2_attenuated) {
  cat("H2 Supported: YES ✓\n")
  cat("Interpretation: Education's negative relationship is substantially reduced\n")
  cat("when cognitive ability (vocabulary) is controlled.\n")
} else if (h2_negative_effect & !h2_attenuated) {
  cat("H2 Supported: PARTIAL\n")
  cat("Interpretation: Negative education effect present but attenuation is modest.\n")
} else {
  cat("H2 Supported: NO\n")
}

write_csv(
  tibble(
    Model = c("Model 1 (demographics)", "Model 2 (+ personality)", "Model 3 (+ cognitive)"),
    Education_beta = round(c(edu_beta_model1, edu_beta_model2, edu_beta_model3), 4),
    Education_b    = round(c(edu_b_model1, edu_b_model2, edu_b_model3), 4),
    Education_p    = c(format_p(edu_p_model1, 4), format_p(edu_p_model2, 4), format_p(edu_p_model3, 4)),
    Attenuation_pct = c(NA,
                        round(((abs(edu_beta_model1) - abs(edu_beta_model2)) /
                                 abs(edu_beta_model1)) * 100, 1),
                        round(attenuation_from_baseline, 1))
  ),
  here("output", "tables", "07b_education_attenuation_analysis.csv")
)
log_step("Education attenuation analysis saved")

# Model 4: Full model (+ Social)
cat("\nModel 4: Full Model (+ Social Variables)\n")
model4 <- lm(build_formula("gcbs_total", preds_m4), data = regression_data)
model4_summary <- summary(model4)
model4_std <- lm.beta(model4)
print_model_fit(model4_summary, model3_summary)
anova_3_4 <- anova(model3, model4)
print_f_change(anova_3_4)

# ---- 7.3 Final Model Detailed Results ----

# FINAL MODEL (Model 4) DETAILED RESULTS

model4_coefs <- cbind(coef(summary(model4)), confint(model4)) %>%
  as.data.frame() %>%
  rownames_to_column("Predictor") %>%
  rename(b = Estimate, SE = `Std. Error`, t = `t value`, p = `Pr(>|t|)`,
         CI_lower = `2.5 %`, CI_upper = `97.5 %`)

betas <- coef(model4_std)[-1]

model4_coefs <- model4_coefs %>%
  filter(Predictor != "(Intercept)") %>%
  mutate(
    beta = betas[match(Predictor, names(betas))],
    sig = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", TRUE ~ ""),
    r = sapply(Predictor, function(pred) {
      if (pred %in% colnames(regression_cor_matrix)) regression_cor_matrix["gcbs_total", pred] else NA_real_
    }),
    beta_minus_r = abs(beta) - abs(r)
  )

print(model4_coefs %>% dplyr::select(Predictor, b, SE, beta, r, t, p, sig))
log_step("Full regression results calculated")

# ---- 7.3.1 H3: Age Linear Negative Relationship ----

# H3 TEST: Age Linear Negative Relationship

age_beta <- coef(model4_std)["age"]
age_b    <- coef(model4)["age"]
age_se   <- model4_summary$coefficients["age", "Std. Error"]
age_t    <- model4_summary$coefficients["age", "t value"]
age_p    <- model4_summary$coefficients["age", "Pr(>|t|)"]
age_r    <- regression_cor_matrix["gcbs_total", "age"]
age_ci   <- confint(model4)["age", ]

cat(sprintf("\nZero-order r: %.4f\n", age_r))
cat(sprintf("b = %.4f, 95%% CI [%.4f, %.4f]\n", age_b, age_ci[1], age_ci[2]))
cat(sprintf("β = %.4f, SE = %.4f, t = %.2f, p = %s\n\n",
            age_beta, age_se, age_t, format_p(age_p)))

h3_negative    <- age_beta < 0
h3_significant <- age_p < 0.05

# H3 RESULT
cat("Negative:", h3_negative, "| Significant:", h3_significant, "\n")
cat("H3 Supported:", ifelse(h3_negative & h3_significant, "YES ✓", "NO"), "\n")

if (h3_negative & h3_significant) {
  cat(sprintf("\nFor each 1-year increase in age, GCBS decreases by %.4f (β = %.3f),\n",
              abs(age_b), age_beta))
  cat("controlling for all other predictors.\n")
}

# Age cohort visualization
regression_data <- regression_data %>%
  mutate(
    age_cohort = cut(age, breaks = c(12, 19, 29, 39, 49, 59, Inf),
                     labels = c("13-19", "20-29", "30-39", "40-49", "50-59", "60+"))
  )

cohort_summary <- regression_data %>%
  group_by(age_cohort) %>%
  summarise(
    n = n(),
    mean_gcbs   = mean(gcbs_total, na.rm = TRUE),
    sd_gcbs     = sd(gcbs_total, na.rm = TRUE),
    se_gcbs     = sd_gcbs / sqrt(n),
    ci_lower    = mean_gcbs - 1.96 * se_gcbs,
    ci_upper    = mean_gcbs + 1.96 * se_gcbs,
    median_gcbs = median(gcbs_total, na.rm = TRUE),
    .groups = "drop"
  )

# GCBS by Age Cohort:
print(cohort_summary)

h3_cohort_plot <- ggplot(cohort_summary, aes(x = age_cohort, y = mean_gcbs)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, linewidth = 0.8) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.5, size = 3, color = "gray30") +
  geom_hline(yintercept = mean(regression_data$gcbs_total, na.rm = TRUE),
             linetype = "dashed", color = "red", alpha = 0.7) +
  labs(title = "Conspiracy Beliefs by Age Cohort",
       subtitle = "H3: Younger cohorts show higher GCBS scores",
       x = "Age Cohort", y = "Mean GCBS Score",
       caption = sprintf("Error bars = 95%% CI | Red dashed = grand mean (%.2f)",
                         mean(regression_data$gcbs_total, na.rm = TRUE))) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40"),
        axis.title = element_text(face = "bold"))

ggsave(here("output", "figures", "14_h3_age_cohort_comparison.png"),
       h3_cohort_plot, width = 10, height = 8, dpi = 300)

write_csv(cohort_summary, here("output", "tables", "07c_age_cohort_summary.csv"))
log_step("Age cohort analysis saved")

# ---- 7.4 Hierarchical Summary Table ----

hier_summary <- tibble(
  Model = c("1. Demographics", "2. + Personality", "3. + Cognitive", "4. + Social"),
  R2       = c(model1_summary$r.squared, model2_summary$r.squared,
               model3_summary$r.squared, model4_summary$r.squared),
  `Adj. R2` = c(model1_summary$adj.r.squared, model2_summary$adj.r.squared,
                model3_summary$adj.r.squared, model4_summary$adj.r.squared),
  `ΔR2` = c(model1_summary$r.squared,
            model2_summary$r.squared - model1_summary$r.squared,
            model3_summary$r.squared - model2_summary$r.squared,
            model4_summary$r.squared - model3_summary$r.squared),
  `F change` = c(model1_summary$fstatistic[1],
                 anova_1_2$F[2], anova_2_3$F[2], anova_3_4$F[2]),
  p = c(pf(model1_summary$fstatistic[1], model1_summary$fstatistic[2],
           model1_summary$fstatistic[3], lower.tail = FALSE),
        anova_1_2$`Pr(>F)`[2], anova_2_3$`Pr(>F)`[2], anova_3_4$`Pr(>F)`[2])
) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

print(hier_summary)
write_csv(hier_summary, here("output", "tables", "06_hierarchical_regression_summary.csv"))
log_step("Hierarchical summary saved")

# ---- 7.5 Suppression Effects Analysis ----

# SUPPRESSION EFFECTS ANALYSIS
# Comparing β to r: Classical suppression when |β| > |r|\n\n

suppressors <- model4_coefs %>%
  filter(
    abs(beta) > abs(r) * SUPPRESSION_EXPLORE |
      (sign(beta) != sign(r) & abs(r) > 0.05)
  ) %>%
  arrange(desc(abs(beta_minus_r)))

if (nrow(suppressors) > 0) {
  cat("SUPPRESSION EFFECTS DETECTED:\n\n")
  print(suppressors %>% dplyr::select(Predictor, r, beta, beta_minus_r, p, sig))
} else {
  log_step("No strong suppression effects detected")
}

write_csv(model4_coefs, here("output", "tables", "07_suppression_analysis.csv"))

# ---- 7.6 Semi-Partial Correlations ----

# SEMI-PARTIAL CORRELATIONS
semipartial_r2 <- effectsize::r2_semipartial(model4)
print(semipartial_r2)

# ---- 7.7 VIF Diagnostics ----

# COLLINEARITY (VIF)

vif_values <- car::vif(model4)
print(round(vif_values, 2))

if (any(vif_values > VIF_SEVERE)) {
  cat("\n SEVERE multicollinearity (VIF > ", VIF_SEVERE, "):\n")
  print(vif_values[vif_values > VIF_SEVERE])
} else if (any(vif_values > VIF_MODERATE)) {
  cat("\n Moderate multicollinearity (VIF > ", VIF_MODERATE, "):\n")
  print(vif_values[vif_values > VIF_MODERATE])
} else {
  log_step(paste("No problematic multicollinearity (all VIF <", VIF_MODERATE, ")"))
}

write_csv(enframe(vif_values, "Predictor", "VIF") %>% arrange(desc(VIF)),
          here("output", "tables", "08_vif_collinearity.csv"))

# ---- 7.8 Save Model Objects ----

saveRDS(model1, here("output", "models", "model1_demographics.rds"))
saveRDS(model2, here("output", "models", "model2_personality.rds"))
saveRDS(model3, here("output", "models", "model3_cognitive.rds"))
saveRDS(model4, here("output", "models", "model4_full.rds"))

log_step("All model objects saved")

# ---- STEP 8: RELATIVE IMPORTANCE ANALYSIS (DOMINANCE ANALYSIS) ----
# LMG metric averages each predictor's R² contribution across all possible
# model orderings — more informative than β when predictors are correlated.

# RELATIVE IMPORTANCE ANALYSIS

# ---- 8.1 LMG Calculation ----

rel_imp <- calc.relimp(model4, type = c("lmg"), rela = TRUE)

rel_imp_df <- data.frame(
  Predictor     = names(rel_imp$lmg),
  LMG           = rel_imp$lmg,
  Percent_of_R2 = (rel_imp$lmg / sum(rel_imp$lmg)) * 100
) %>%
  arrange(desc(LMG))

# Relative Importance (LMG):
print(rel_imp_df)
cat("\nTotal R² (model):", round(model4_summary$r.squared, 4), "\n")
cat("LMG proportions sum to:", round(sum(rel_imp$lmg), 4), "(normalized)\n")

# ---- 8.2 Visualization ----

rel_imp_plot <- ggplot(rel_imp_df, aes(x = reorder(Predictor, LMG), y = LMG)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = paste0(round(Percent_of_R2, 1), "%")),
            hjust = -0.2, size = 3) +
  coord_flip() +
  labs(title = "Relative Importance of Predictors (LMG Metric)",
       subtitle = paste0("Total R² = ", round(model4_summary$r.squared, 4)),
       x = "Predictor", y = "Proportion of R²",
       caption = "Averaged across all possible model orderings") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 10))

ggsave(here("output", "figures", "07_relative_importance.png"),
       rel_imp_plot, width = 10, height = 8, dpi = 300)
log_step("Relative importance plot saved")

# ---- 8.3 Compare LMG vs Standardized Betas ----

comparison <- model4_coefs %>%
  dplyr::select(Predictor, beta, p) %>%
  left_join(rel_imp_df, by = "Predictor") %>%
  mutate(abs_beta = abs(beta), rank_beta = rank(-abs_beta), rank_lmg = rank(-LMG)) %>%
  arrange(desc(LMG))

cat("\n--- β vs LMG Rankings ---\n")
print(comparison %>% dplyr::select(Predictor, beta, LMG, Percent_of_R2, rank_beta, rank_lmg))

write_csv(comparison, here("output", "tables", "09_relative_importance_comparison.csv"))
log_step("Comparison table saved")

# ---- STEP 9: REGRESSION ASSUMPTIONS & DIAGNOSTICS ----

# REGRESSION ASSUMPTIONS

# ---- 9.1 Linearity ----

# Linearity (Component-Residual Plots)

png(here("output", "figures", "08_component_residual_plots.png"),
    width = 12, height = 10, units = "in", res = 300)
car::crPlots(model4)
dev.off()

log_step("CR plots saved (check pink vs blue line deviation)")

# ---- 9.2 Independence ----

# Independence (Durbin-Watson)

dw_test <- car::durbinWatsonTest(model4)
cat("DW =", round(dw_test$dw, 3), "| p =", format_p(dw_test$p), "\n")
cat(ifelse(dw_test$dw > 1.5 & dw_test$dw < 2.5,
           "Independence satisfied (DW ≈ 2)\n",
           "Potential autocorrelation\n"))

# ---- 9.3 Homoscedasticity ----

# Homoscedasticity (Breusch-Pagan)

bp_test <- lmtest::bptest(model4)
cat("BP =", round(bp_test$statistic, 3), "| p =", format_p(bp_test$p.value), "\n")

if (bp_test$p.value < 0.05) {
  cat("Heteroscedasticity detected — see robust SEs in Step 10\n")
} else {
  log_step("Homoscedasticity satisfied")
}

resid_fitted_plot <- ggplot(data = NULL, aes(x = fitted(model4), y = residuals(model4))) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(se = TRUE, color = "blue") +
  labs(title = "Residuals vs Fitted", subtitle = "Checking homoscedasticity",
       x = "Fitted Values", y = "Residuals",
       caption = "Ideal: random scatter with constant spread") +
  theme_minimal()

ggsave(here("output", "figures", "09_residuals_vs_fitted.png"),
       resid_fitted_plot, width = 10, height = 8, dpi = 300)

# ---- 9.4 Normality of Residuals ----

# Shapiro-Wilk (subsample if n is large — test is very sensitive)
if (n_regression > SHAPIRO_SUBSAMPLE_N) {
  set.seed(SEED)
  resid_sample <- sample(residuals(model4), SHAPIRO_SUBSAMPLE_N)
  sw_test <- shapiro.test(resid_sample)
  cat("(Subsampled", SHAPIRO_SUBSAMPLE_N, "residuals for Shapiro-Wilk)\n")
} else {
  sw_test <- shapiro.test(residuals(model4))
}

cat("W =", round(sw_test$statistic, 4), "| p =", format_p(sw_test$p.value), "\n")
cat("Note: With large n, SW is very sensitive; visual inspection is more informative.\n")

qq_plot <- ggplot(data = NULL, aes(sample = residuals(model4))) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot: Residuals", x = "Theoretical Quantiles",
       y = "Sample Quantiles") + theme_minimal()

resid_df <- data.frame(resid = residuals(model4))

resid_hist <- ggplot(resid_df, aes(x = resid)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "skyblue", alpha = 0.7, color = "white") +
  geom_density(color = "darkblue", linewidth = 1) +
  stat_function(fun = dnorm,
                args = list(mean = mean(resid_df$resid), sd = sd(resid_df$resid)),
                color = "red", linetype = "dashed", linewidth = 1) +
  labs(title = "Residual Distribution", x = "Residuals", y = "Density") +
  theme_minimal()

normality_combined <- qq_plot / resid_hist +
  plot_annotation(title = "Normality of Residuals",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(here("output", "figures", "10_normality_diagnostics.png"),
       normality_combined, width = 10, height = 10, dpi = 300)
log_step("Normality plots saved")

# ---- 9.5 Influential Cases ----

# Influential Cases & Outliers

influence_stats <- data.frame(
  case      = seq_len(n_regression),
  cooks_d   = cooks.distance(model4),
  dffits    = dffits(model4),
  dfbetas_max = apply(abs(dfbetas(model4)), 1, max),
  leverage  = hatvalues(model4),
  std_resid = rstandard(model4)
)

k <- length(coef(model4)) - 1
cooks_cutoff   <- COOKS_D_MULTIPLIER / n_regression
dffits_cutoff  <- 2 * sqrt(k / n_regression)
leverage_cutoff <- 3 * (k + 1) / n_regression

n_influential_cooks  <- sum(influence_stats$cooks_d > cooks_cutoff)
n_influential_dffits <- sum(abs(influence_stats$dffits) > dffits_cutoff)
n_high_leverage      <- sum(influence_stats$leverage > leverage_cutoff)
n_extreme_resid      <- sum(abs(influence_stats$std_resid) > 3)

cat(sprintf("Cook's D > %.5f: %d (%.2f%%)\n",
            cooks_cutoff, n_influential_cooks, n_influential_cooks / n_regression * 100))
cat(sprintf("DFFITS > %.3f: %d (%.2f%%)\n",
            dffits_cutoff, n_influential_dffits, n_influential_dffits / n_regression * 100))
cat(sprintf("Leverage > %.4f: %d (%.2f%%)\n",
            leverage_cutoff, n_high_leverage, n_high_leverage / n_regression * 100))
cat(sprintf("|Std residual| > 3: %d (%.2f%%)\n",
            n_extreme_resid, n_extreme_resid / n_regression * 100))

influence_plot <- ggplot(influence_stats, aes(x = leverage, y = std_resid)) +
  geom_point(aes(size = cooks_d, color = cooks_d > cooks_cutoff), alpha = 0.5) +
  geom_hline(yintercept = c(-3, 3), color = "red", linetype = "dashed") +
  geom_vline(xintercept = leverage_cutoff, color = "blue", linetype = "dashed") +
  scale_color_manual(values = c("grey40", "red"),
                     labels = c("Not influential", "Influential")) +
  labs(title = "Influence Plot", x = "Leverage", y = "Standardized Residuals",
       size = "Cook's D", color = "Status",
       caption = sprintf("Red dashed = ±3 SD | Blue dashed = leverage (%.4f)",
                         leverage_cutoff)) +
  theme_minimal() + theme(legend.position = "right")

ggsave(here("output", "figures", "11_influence_plot.png"),
       influence_plot, width = 10, height = 8, dpi = 300)
log_step("Influence plot saved")

# ---- 9.6 Assumption Summary ----

assumption_summary <- tibble(
  Assumption = c("Linearity", "Independence", "Homoscedasticity",
                 "Normality", "No Influential Outliers"),
  Status = c(
    "Check plots",
    ifelse(dw_test$dw > 1.5 & dw_test$dw < 2.5, "✓ Satisfied", "⚠ Concern"),
    ifelse(bp_test$p.value >= 0.05, "✓ Satisfied", "⚠ Violated"),
    "Check plots",
    paste0(n_influential_cooks, " influential cases")
  ),
  Test = c("—", paste0("DW = ", round(dw_test$dw, 3)),
           paste0("BP = ", round(bp_test$statistic, 3)),
           paste0("W = ", round(sw_test$statistic, 4)),
           paste0("Cook's D > ", round(cooks_cutoff, 5))),
  P = c("—", format_p(dw_test$p), format_p(bp_test$p.value),
        format_p(sw_test$p.value), "—")
)

print(assumption_summary)
write_csv(assumption_summary, here("output", "tables", "10_assumption_check_summary.csv"))
log_step("Assumption summary saved")

# ---- STEP 10: ROBUSTNESS & SENSITIVITY ANALYSES ----

# ROBUSTNESS CHECKS

# ---- 10.1 Exclude Influential Cases ----

# Sensitivity: Excluding Influential Cases

influential_cases <- which(influence_stats$cooks_d > cooks_cutoff)

if (length(influential_cases) > 0 & length(influential_cases) < nrow(regression_data) * 0.05) {
  cat("Removing", length(influential_cases), "cases with Cook's D >",
      round(cooks_cutoff, 5), "\n\n")
  
  regression_data_robust <- regression_data[-influential_cases, ]
  model4_robust <- lm(build_formula("gcbs_total", preds_m4), data = regression_data_robust)
  model4_robust_summary <- summary(model4_robust)
  
  cat("Full sample:   N =", n_regression,
      "| R² =", round(model4_summary$r.squared, 4), "\n")
  cat("Robust sample: N =", nrow(regression_data_robust),
      "| R² =", round(model4_robust_summary$r.squared, 4), "\n\n")
  
  coef_comparison <- data.frame(
    Predictor    = names(coef(model4))[-1],
    Beta_Full    = coef(lm.beta(model4))[-1],
    Beta_Robust  = coef(lm.beta(model4_robust))[-1]
  ) %>%
    mutate(Difference = Beta_Robust - Beta_Full,
           Pct_Change = ifelse(abs(Beta_Full) < 0.001, NA_real_,
                               (Difference / Beta_Full) * 100)) %>%
    arrange(desc(abs(Pct_Change)))
  
  cat("Largest coefficient changes:\n")
  print(head(coef_comparison, 5))
  
  write_csv(coef_comparison, here("output", "tables", "11_robustness_coefficient_comparison.csv"))
  
  cat("\nResults appear",
      ifelse(max(abs(coef_comparison$Pct_Change), na.rm = TRUE) < 20,
             "robust (changes < 20%)", "somewhat sensitive to influential cases"), "\n")
  
} else if (length(influential_cases) == 0) {
  log_step("No influential cases — robust analysis not needed")
} else {
  cat("Too many influential cases (>5%) — examine data quality\n")
}

# ---- 10.2 Robust Standard Errors ----

# Robust Standard Errors (HC3)

robust_se <- lmtest::coeftest(model4, vcov = sandwich::vcovHC(model4, type = "HC3"))
print(robust_se)

robust_results <- as.data.frame(unclass(robust_se), check.names = FALSE) %>%
  rownames_to_column("Predictor") %>%
  rename(b = Estimate, SE_robust = `Std. Error`, t_robust = `t value`, p_robust = `Pr(>|t|)`)

write_csv(robust_results, here("output", "tables", "12_robust_standard_errors.csv"))
log_step("Robust SEs saved (use if BP test indicated heteroscedasticity)")

# ---- 10.3 Subgroup: Gender ----

# Subgroup: Gender

data_male   <- regression_data %>% filter(gender == 0)
data_female <- regression_data %>% filter(gender == 1)

# Subgroup formula omits gender from predictors
preds_subgroup <- setdiff(preds_m4, "gender")

if (nrow(data_male) >= 100 & nrow(data_female) >= 100) {
  model_male   <- lm(build_formula("gcbs_total", preds_subgroup), data = data_male)
  model_female <- lm(build_formula("gcbs_total", preds_subgroup), data = data_female)
  
  cat(sprintf("Males:   N = %d, R² = %.4f\n", nrow(data_male), summary(model_male)$r.squared))
  cat(sprintf("Females: N = %d, R² = %.4f\n\n", nrow(data_female), summary(model_female)$r.squared))
  
  for (pred in c("education_num", "openness", "agreeableness")) {
    if (pred %in% names(coef(model_male))) {
      cat(sprintf("  %s: β_male = %.3f | β_female = %.3f\n",
                  pred, coef(lm.beta(model_male))[pred], coef(lm.beta(model_female))[pred]))
    }
  }
} else {
  cat("Insufficient N in one or both groups for subgroup analysis\n")
}

# ---- 10.4 Subgroup: Age Groups ----

# Subgroup: Age Groups

regression_data <- regression_data %>%
  mutate(age_group = case_when(
    age < 25 ~ "Young (< 25)",
    age < 40 ~ "Middle (25–39)",
    TRUE      ~ "Older (40+)"
  ))

print(regression_data %>%
        group_by(age_group) %>%
        summarise(n = n(), mean_gcbs = mean(gcbs_total), sd_gcbs = sd(gcbs_total),
                  .groups = "drop"))

# ---- COMPREHENSIVE HYPOTHESIS TESTING SUMMARY ----

# COMPREHENSIVE HYPOTHESIS TESTING SUMMARY

hypothesis_results <- tibble(
  Hypothesis = c(
    "H1a: Agreeableness strongest zero-order personality predictor",
    "H1b: Openness positive suppression after cognitive control",
    "H2: Education effect attenuated by vocabulary",
    "H3: Age demonstrates significant negative linear relationship"
  ),
  Prediction = c(
    "Agreeableness has largest |r| with GCBS among Big Five",
    "Openness β becomes positive/larger after vocab (|β| > |r|)",
    "Education β reduces >20% when vocabulary added",
    "Age β < 0 and p < .05 in final model"
  ),
  Key_Statistic = c(
    sprintf("Agreeableness r = %.3f (rank %d of 5)", agreeableness_r, agreeableness_rank),
    sprintf("Openness: r = %.3f → β = %.3f", openness_r, openness_beta_model3),
    sprintf("Education β attenuation = %.1f%%", attenuation_pct),
    sprintf("Age β = %.3f, p = %s", age_beta, format_p(age_p))
  ),
  Supported = c(
    ifelse(h1a_supported, "Yes ✓", "No"),
    ifelse(suppression_detected & positive_after_control, "Yes ✓", "No"),
    ifelse(h2_negative_effect & h2_attenuated, "Yes ✓",
           ifelse(h2_negative_effect, "Partial", "No")),
    ifelse(h3_negative & h3_significant, "Yes ✓", "No")
  )
)

# HYPOTHESIS TESTING RESULTS:
print(hypothesis_results, n = Inf, width = Inf)

n_supported     <- sum(hypothesis_results$Supported == "Yes ✓")
n_partial       <- sum(hypothesis_results$Supported == "Partial")
n_not_supported <- sum(hypothesis_results$Supported == "No")

cat(sprintf("\nSUMMARY: %d fully supported, %d partial, %d not supported\n",
            n_supported, n_partial, n_not_supported))

write_csv(hypothesis_results, here("output", "tables", "14_hypothesis_testing_summary.csv"))
log_step("Hypothesis testing summary saved")

# ---- STEP 11: PUBLICATION-READY OUTPUTS ----

# PUBLICATION-READY OUTPUTS

# ---- 11.1 APA Regression Table ----

pub_table <- model4_coefs %>%
  mutate(Predictor_clean = predictor_labels[Predictor]) %>%
  dplyr::select(Predictor_clean, b, SE, beta, t, p) %>%
  mutate(
    b    = sprintf("%.3f", b),
    SE   = sprintf("%.3f", SE),
    beta = sprintf("%.3f", beta),
    t    = sprintf("%.2f", t),
    p    = ifelse(p < .001, "<.001", sprintf("%.3f", p))
  ) %>%
  rename(Predictor = Predictor_clean, B = b, `β` = beta)

write_csv(pub_table, here("output", "tables", "13_regression_table_APA.csv"))

kable(pub_table,
      caption = "Table 2. Hierarchical Multiple Regression Predicting Conspiracy Beliefs") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  footnote(general = sprintf("N = %d. R² = %.3f, Adj. R² = %.3f, F(%g, %g) = %.2f, p < .001",
                             n_regression, model4_summary$r.squared,
                             model4_summary$adj.r.squared,
                             model4_summary$fstatistic[2], model4_summary$fstatistic[3],
                             model4_summary$fstatistic[1])) %>%
  save_kable(here("output", "tables", "13_regression_table_APA.html"))

log_step("APA regression table saved")

# ---- 11.2 Coefficient Plot ----

std_params <- effectsize::standardize_parameters(model4, method = "refit")

std_coef_plot_data <- std_params %>%
  filter(Parameter != "(Intercept)") %>%
  left_join(model4_coefs %>% dplyr::select(Predictor, p), 
            by = c("Parameter" = "Predictor")) %>%
  mutate(
    Predictor_clean = predictor_labels[Parameter],
    significant = p < 0.05
  )

coef_plot <- std_coef_plot_data %>%
  ggplot(aes(x = reorder(Predictor_clean, Std_Coefficient),
             y = Std_Coefficient, color = significant)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  scale_color_manual(values = c("grey60", "darkblue"), labels = c("ns", "p < .05")) +
  labs(title = "Standardized Coefficients (β) with 95% CIs",
       subtitle = "Predicting Conspiracy Beliefs (GCBS)",
       x = "Predictor", y = "Standardized β", color = "Significance") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")

ggsave(here("output", "figures", "12_coefficient_plot.png"),
       coef_plot, width = 10, height = 8, dpi = 300)
log_step("Coefficient plot saved")

# ---- 11.3 Model Comparison Plot ----

model_comparison_plot <- hier_summary %>%
  ggplot(aes(x = Model, y = R2, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(size = 4, color = "darkblue") +
  geom_text(aes(label = sprintf("R² = %.3f", R2)), vjust = -1, size = 3.5) +
  ylim(0, max(hier_summary$R2) * 1.15) +
  labs(title = "Hierarchical Regression: Incremental R²",
       subtitle = "Each block adds to prediction of conspiracy beliefs",
       x = "Model", y = "R²") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here("output", "figures", "13_model_comparison.png"),
       model_comparison_plot, width = 10, height = 6, dpi = 300)
log_step("Model comparison plot saved")

# ---- 11.4 Results Section Template ----

results_template <- paste0(
  "RESULTS SECTION TEMPLATE\n", strrep("=", 70), "\n\n",
  
  "Sample Characteristics\n", strrep("-", 70), "\n",
  "After validity screening and listwise deletion, N = ", n_regression, " participants (",
  round(sum(regression_data$gender == 1) / n_regression * 100, 1),
  "% female). Mean age = ", round(mean(regression_data$age), 1),
  " (SD = ", round(sd(regression_data$age), 1), ").\n\n",
  
  "Reliability\n", strrep("-", 70), "\n",
  "GCBS: excellent internal consistency (α = ", round(alpha_result$total$std.alpha, 2),
  ", ω = ", round(omega_result$omega.tot, 2), ").\n\n",
  
  "Descriptive Statistics\n", strrep("-", 70), "\n",
  "GCBS: ", round(min(gcbs_analysis$gcbs_total, na.rm = TRUE), 2), "–",
  round(max(gcbs_analysis$gcbs_total, na.rm = TRUE), 2),
  " (M = ", round(mean(gcbs_analysis$gcbs_total, na.rm = TRUE), 2),
  ", SD = ", round(sd(gcbs_analysis$gcbs_total, na.rm = TRUE), 2), ").\n\n",
  
  "Hierarchical Regression\n", strrep("-", 70), "\n",
  "Final model: F(", model4_summary$fstatistic[2], ", ",
  model4_summary$fstatistic[3], ") = ", round(model4_summary$fstatistic[1], 2),
  ", p < .001, R² = ", round(model4_summary$r.squared, 3),
  ", Adj. R² = ", round(model4_summary$adj.r.squared, 3), ".\n\n",
  
  "Hypothesis Testing & Theoretical Integration\n", strrep("-", 70), "\n\n",
  "H1a (Interpersonal Layer): ", strongest_personality, " emerged as the strongest zero-order predictor (r = ", 
  round(personality_cors$r[1], 3), "). ", 
  ifelse(h1a_supported, 
         "H1a supported. This aligns with the premise that trait antagonism and epistemic mistrust form the baseline vulnerability for conspiratorial ideation.",
         paste0("H1a not supported; Agreeableness ranked ", agreeableness_rank, " instead.")), "\n\n",
  
  "H1b & H2 (Cognitive Layer - Dual-Process Theory): \n",
  "Openness ",
  ifelse(suppression_detected & positive_after_control,
         paste0("showed a classic suppression effect, with its β increasing from ", round(openness_beta_model2, 3),
                " to ", round(openness_beta_model3, 3), " after controlling for cognitive ability. H1b supported."),
         paste0("did not show clear suppression. H1b not supported.")), "\n",
  "Furthermore, the Education β ",
  ifelse(h2_negative_effect, paste0("was initially negative (Model 1 β = ", round(edu_beta_model1, 3), "). "), "was not significant. "),
  ifelse(h2_attenuated,
         paste0("It attenuated by ", round(attenuation_pct, 1), "% when vocabulary was added. H2 supported."),
         paste0("The attenuation of ", round(attenuation_pct, 1), "% did not meet the threshold. H2 ",
                ifelse(h2_negative_effect, "partially", "not"), " supported.")), "\n",
  ifelse((suppression_detected & positive_after_control) | (h2_negative_effect & h2_attenuated),
         "Together, these findings support a dual-process framework: raw cognitive capacity (System 2) acts as a stronger epistemic buffer against apophenia than formal education, while unchecked Openness (System 1) increases conspiratorial pattern-seeking.\n\n",
         "The dual-process framework was not clearly supported by the current data; neither the suppression nor attenuation criteria were fully met.\n\n"),
  
  "H3 (Environmental Layer - Digital Cohorts): Age demonstrated a significant linear effect (β = ", 
  round(age_beta, 3), ", p = ", format_p(age_p), "). ",
  ifelse(h3_negative & h3_significant, 
         "H3 supported. This bolsters the digital cohort effect model, suggesting younger generations immersed in decentralized digital environments exhibit higher baseline vulnerability to these narratives.", 
         "H3 not supported."), "\n\n",
  
  "Assumption Checks\n", strrep("-", 70), "\n",
  "[INSERT: Summary from diagnostics section]\n\n",
  
  "Robustness\n", strrep("-", 70), "\n",
  "[INSERT: Summary from robustness analyses]\n\n"
)

writeLines(results_template, here("documentation", "results_section_template.txt"))
log_step("Results template saved with integrated theoretical frameworks")

# ---- STEP 12: FINAL DOCUMENTATION & REPRODUCIBILITY ----

# FINAL DOCUMENTATION

# ---- 12.1 Analysis Summary ----

analysis_summary <- list(
  dataset_info = list(
    name = "GCBS (Generic Conspiracist Beliefs Scale)",
    source = "Open Psychometrics, 2016",
    raw_n = nrow(gcbs_raw),
    final_n = n_final,
    exclusion_rate = round((1 - n_final / nrow(gcbs_raw)) * 100, 2)
  ),
  reliability = list(
    cronbach_alpha = round(alpha_result$total$std.alpha, 3),
    mcdonald_omega = round(omega_result$omega.tot, 3)
  ),
  regression_results = list(
    final_r_squared = round(model4_summary$r.squared, 4),
    adj_r_squared   = round(model4_summary$adj.r.squared, 4),
    f_statistic     = round(model4_summary$fstatistic[1], 2),
    p_value         = format_p(pf(model4_summary$fstatistic[1],
                                  model4_summary$fstatistic[2],
                                  model4_summary$fstatistic[3],
                                  lower.tail = FALSE))
  ),
  top_predictors = model4_coefs %>%
    arrange(p) %>% head(5) %>%
    dplyr::select(Predictor, beta, p) %>% as.list()
)

saveRDS(analysis_summary, here("output", "analysis_summary.rds"))
log_step("Analysis summary saved")

# ---- 12.2 Session Info (End) ----

sink(here("documentation", "session_info_end.txt"))
print(sessionInfo())
sink()
log_step("End session info saved")

# ---- 12.3 Completion Log ----

completion_log <- tibble(
  timestamp       = Sys.time(),
  script          = "GCBS_hierarchical_regression_analysis.R",
  status          = "Completed successfully",
  sample_size     = n_final,
  r_version       = paste(R.version$major, R.version$minor, sep = "."),
  models_saved    = 4,
  tables_created  = 17,
  figures_created = 14
)

write_csv(completion_log, here("documentation", "analysis_completion_log.csv"))
log_step("Completion log saved")

# ---- ANALYSIS COMPLETE ---- 

# ANALYSIS COMPLETE

# Summary:
cat("  Cleaned", n_final, "cases from", nrow(gcbs_raw), "raw observations\n")
cat("  GCBS reliability: α =", round(alpha_result$total$std.alpha, 3),
    "| ω =", round(omega_result$omega.tot, 3), "\n")
cat("  Final model R² =", round(model4_summary$r.squared, 4), "\n")
cat("  Hypotheses:", n_supported, "supported,", n_partial, "partial,",
    n_not_supported, "not supported\n")
cat("  14 figures, 17 tables, 4 model objects saved\n\n")

cat("Output locations:\n")
cat("  Tables:        ", here("output", "tables"), "\n")
cat("  Figures:       ", here("output", "figures"), "\n")
cat("  Models:        ", here("output", "models"), "\n")
cat("  Documentation: ", here("documentation"), "\n\n")

# Next steps:
# 1. Review diagnostic plots
# 2. Write up results using the template
# 3. Interpret suppression effects theoretically
# 4. Prepare manuscript tables and figures

# End of Script

# ---- END OF GCBS CONSPIRACY BELIEFS ANALYSIS ----