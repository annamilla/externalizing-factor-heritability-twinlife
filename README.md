# Externalizing Factor and Heritability Across Age and Gender in TwinLife

Analysis code for deriving a general externalizing factor from TwinLife data and estimating its heritability with twin models.

## Contents

- `R/` — data preparation, factor analysis, invariance tests, age/gender analyses, and Mplus export
- `models/` — lavaan latent externalizing factor models and sensitivity analyses
- `Mplus/` — ACE, ADE, and AE twin-model input files
- `onyx/` — path diagram of latent externalizing factor structure 

## Requirements

- R
- R packages: `dplyr`, `tidyr`, `haven`, `ggplot2`, `Hmisc`, `lavaan`, `semTools`, `semPlot`, `corrplot`
- Mplus for the twin models

## Workflow

1. Replace the local paths defined near the top of the scripts.
2. Run the scripts in R in numerical order.
3. Run .inp files in Mplus.


