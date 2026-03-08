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
RAPID_RESPONSE_SD       <- 2       # Flag respondents faster than M - 2*SD
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