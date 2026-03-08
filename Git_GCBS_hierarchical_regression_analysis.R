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

