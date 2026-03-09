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
#
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
#
# SAMPLE SIZE & POWER (Field et al., 2012, p. 273):
#   With ~15 predictors: minimum N = 50 + 8(15) = 170
#   For medium effect (f² = .15), power = .80: ~150 cases
#   Expected dataset: several thousand (more than adequate)

# STEP 1: PROJECT SETUP & CONFIGURATION

# 1.1 Analysis Parameters
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

# STEP 2: DATA IMPORT & INITIAL EXPLORATION

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
# STEP 3: DATA CLEANING & PREPARATION

# ---- 3.1 Validity Screening ----

# 3.1.1 Vocabulary check: flag participants endorsing ≥2 fake words (VCL6, VCL9, VCL12)
gcbs_clean <- gcbs_raw %>%
  mutate(
    fake_words_checked = rowSums(pick(VCL6, VCL9, VCL12), na.rm = TRUE),
    invalid_vocabulary = fake_words_checked >= FAKE_WORD_THRESHOLD
  )

n_invalid_vocab <- sum(gcbs_clean$invalid_vocabulary, na.rm = TRUE)
cat("\n=== Vocabulary Validity ===\n")
cat("Cases endorsing ≥2 fake words:", n_invalid_vocab,
    sprintf("(%.2f%%)\n", n_invalid_vocab / nrow(gcbs_clean) * 100))

# 3.1.2 Age validity: exclude implausible ages or missing
gcbs_clean <- gcbs_clean %>%
  mutate(invalid_age = age < AGE_RANGE[1] | age > AGE_RANGE[2] | is.na(age))

cat("\n=== Age Validity ===\n")
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

cat("\n=== Response Time Validity ===\n")
cat("Median:", round(time_median, 1), "s | IQR:", round(time_iqr, 1), "s\n")
cat("Calculated Bound (Med - 1.5*IQR):", round(time_cutoff, 1), "s\n")
cat("Applied Cutoff (Max of bound or 30s):", round(final_cutoff, 1), "s\n")
cat("Flagged as too fast:", sum(gcbs_clean$rapid_responding, na.rm = TRUE), "\n")

# 3.1.4 Combine all validity flags
gcbs_clean <- gcbs_clean %>%
  mutate(valid_case = !(invalid_vocabulary | invalid_age | rapid_responding))

n_excluded <- sum(!gcbs_clean$valid_case, na.rm = TRUE)
cat("\n=== Overall Validity ===\n")
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

cat("\n=== GCBS Scoring ===\n")
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

cat("\nTIPI Summary:\n")
print(tipi_summary)

# ---- 3.4 Additional Predictor Variables ----

# 3.4.1 Vocabulary score (real words only, excluding fakes VCL6/9/12)
real_words <- paste0("VCL", c(1:5, 7:8, 10:11, 13:16))

gcbs_clean <- gcbs_clean %>%
  mutate(
    vocab_score       = rowSums(pick(all_of(real_words)), na.rm = TRUE),
    vocab_total_items = rowSums(!is.na(pick(all_of(real_words))))
  )

cat("\n=== Vocabulary Score ===\n")
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
