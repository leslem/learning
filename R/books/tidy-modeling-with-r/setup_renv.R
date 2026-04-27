renv::settings$snapshot.type("all")
options(renv.config.pak.enabled = TRUE)

renv::init(bare = TRUE)

my_repos <- c(
  'Posit CRAN' = "https://packagemanager.posit.co/cran/latest"
)
options(repos = my_repos)
options()$repos

renv::install(
  c(
    "applicable",
    "av",
    "baguette",
    "beans",
    "bestNormalize",
    "bookdown",
    "broom",
    "censored",
    "corrplot",
    "corrr",
    "Cubist",
    "DALEXtra",
    "dials",
    "dimRed",
    "discrim",
    "doMC",
    "dplyr",
    "earth",
    "embed",
    "fastICA",
    "finetune",
    "forcats",
    "ggforce",
    "ggplot2",
    "glmnet",
    "gridExtra",
    "infer",
    "kableExtra",
    "kernlab",
    "kknn",
    "klaR",
    "knitr",
    # "learntidymodels",
    "lime",
    "lme4",
    "lubridate",
    "mda",
    # "mixOmics",
    "modeldata",
    "multilevelmod",
    "nlme",
    "nnet",
    "parsnip",
    "patchwork",
    "pillar",
    "poissonreg",
    "prettyunits",
    "probably",
    "pscl",
    "purrr",
    "ranger",
    "recipes",
    "rlang",
    "rmarkdown",
    "rpart",
    "rsample",
    "rstanarm",
    "rules",
    "sessioninfo",
    "stacks",
    "stringr",
    "svglite",
    "text2vec",
    "textrecipes",
    "themis",
    "tibble",
    "tidymodels",
    "tidyposterior",
    "tidyverse",
    "tune",
    "uwot",
    "workflows",
    "workflowsets",
    "xgboost",
    "yardstick"
  ),
  prompt = F
)

renv::snapshot(type = "all", prompt = FALSE)

renv::restore(prompt = FALSE)
