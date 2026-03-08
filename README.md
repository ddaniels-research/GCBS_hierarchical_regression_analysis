# Psychological and Demographic Antecedents of Conspiracy Beliefs
**A Hierarchical Regression Analysis**

> **WORK IN PROGRESS**
> *This repository is currently under active development.*

## Project Overview
This project examines the personality, cognitive, and demographic predictors of conspiracy beliefs using a hierarchical regression framework. The analysis approaches conspiracy ideation through a unified, multi-layered systems framework:

1. **Interpersonal Layer:** Examining trait antagonism and epistemic mistrust (Agreeableness).
2. **Cognitive Layer:** Exploring Dual-Process Theory and apophenia, specifically how cognitive ability (Vocabulary) might suppress or attenuate the effects of System 1 heuristic exploration (Openness).
3. **Environmental Layer:** Investigating digital cohort effects by examining generational vulnerability to conspiratorial narratives.

**Research Question:** Which personality, cognitive, and demographic factors uniquely predict conspiracy beliefs, and are there suppression effects among these variables?

## Data Source
The data for this analysis comes from the **Generic Conspiracist Beliefs Scale (GCBS)** dataset, collected via online administration in 2016 by Open Psychometrics. 
* *Note: The dataset contains no personally identifiable information, and all participants opted in for research use.*

## Methodology
The core of this project relies on a 4-block Hierarchical Multiple Regression, structured theoretically:
* **Block 1 (Demographics):** Age, gender, education, urban/rural background.
* **Block 2 (Personality):** Big Five dimensions (TIPI).
* **Block 3 (Cognitive):** Vocabulary knowledge (proxy for cognitive ability/System 2 processing).
* **Block 4 (Social):** Marital status, family size.

**Key Analytical Features in the Script:**
* Comprehensive psychometric validation (Cronbach's $\alpha$, McDonald's $\omega$).
* Classical suppression effect testing.
* Relative importance analysis (Dominance Analysis via the LMG metric).
* Robust assumption checking and sensitivity analyses (Cook's D, DFFITS, HC3 robust standard errors).
