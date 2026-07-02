# Code from previous stepslibrary(tidymodels)
library(tidymodels)
data(ames)
ames <- mutate(ames, Sale_Price = log10(Sale_Price))

set.seed(502)
ames_split <- initial_split(ames, prop = 0.80, strata = Sale_Price)
ames_train <- training(ames_split)
ames_test <- testing(ames_split)

ames_rec <-
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + Latitude + Longitude, data = ames_train) |>
  step_log(Gr_Liv_Area, base = 10) |>
  step_other(Neighborhood, threshold = 0.01) |>
  step_dummy(all_nominal_predictors()) |>
  step_interact(~ Gr_Liv_Area:starts_with("Bldg_Type_")) |>
  step_ns(Latitude, Longitude, deg_free = 20)

lm_model <- linear_reg() |> set_engine("lm")

lm_wflow <-
  workflow() |>
  add_model(lm_model) |>
  add_recipe(ames_rec)

lm_fit <- fit(lm_wflow, ames_train)

# Code to fit a random forest model
# no pre-processing required for this model, so no recipe used here
rf_model <-
  rand_forest(trees = 1000) |>
  set_engine("ranger") |>
  set_mode("regression")

rf_wflow <-
  workflow() |>
  add_formula(
    Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + Latitude + Longitude
  ) |>
  add_model(rf_model)

rf_fit <- rf_wflow |> fit(data = ames_train)

# Predict training set to get "Apparent/Resubstitution metric"
# RF is much better on training set, but worse on test set!
estimate_perf <- function(model, dat) {
  # Capture the names of the `model` and `dat` objects
  cl <- match.call()
  obj_name <- as.character(cl$model)
  data_name <- as.character(cl$dat)
  data_name <- gsub("ames_", "", data_name)

  # Estimate these metrics:
  reg_metrics <- metric_set(rmse, rsq)

  model |>
    predict(dat) |>
    bind_cols(dat |> select(Sale_Price)) |>
    reg_metrics(Sale_Price, .pred) |>
    select(-.estimator) |>
    mutate(object = obj_name, data = data_name)
}

estimate_perf(rf_fit, ames_train)
estimate_perf(lm_fit, ames_train)

estimate_perf(rf_fit, ames_test)

# V-fold cross validation
set.seed(1001)
ames_folds <- vfold_cv(ames_train, v = 10)
ames_folds
class(ames_folds$splits[1][[1]])

# You can get to the analysis or assessment sets within a single fold
ames_folds$splits[[1]] |> analysis()
ames_folds$splits[[1]] |> assessment()
# But you generally don't need to since the tidymodels functions will just work with this object

# Do repeats of v-fold CV
ames_r_folds <- vfold_cv(ames_train, v = 10, repeats = 5)
ames_r_folds |> print(n = 100)

# Monte carlo CV
ames_mc_folds <- mc_cv(ames_train, prop = 9 / 10, times = 20)
ames_mc_folds

# Validation set
set.seed(52)
# To put 60% into training, 20% in validation, and 20% in testing:
ames_val_split <- initial_validation_split(ames, prop = c(0.6, 0.2))
ames_val_split

# Object used for resampling:
val_set <- validation_set(ames_val_split)
val_set

# Bootstrapping
ames_boots <- bootstraps(ames_train, times = 5)

# Estimating performance
keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)

set.seed(1003)
rf_res <-
  rf_wflow |>
  fit_resamples(resamples = ames_folds, control = keep_pred)
rf_res
rf_metrics <- collect_metrics(rf_res)
# all replicates
collect_metrics(rf_res, summarize = FALSE)

rf_metrics
# Compare this to the results from the beginning! (64 - 67)
estimate_perf(rf_fit, ames_train)
estimate_perf(lm_fit, ames_train)
estimate_perf(rf_fit, ames_test)

assess_res <- collect_predictions(rf_res)
assess_res
assess_res |>
  ggplot(aes(x = Sale_Price, y = .pred)) +
  geom_point(alpha = .15) +
  geom_abline(color = "red") +
  coord_obs_pred() +
  ylab("Predicted")

# Look into two data points with over-estimates for the predicted sale price
over_predicted <-
  assess_res |>
  mutate(residual = Sale_Price - .pred) |>
  arrange(desc(abs(residual))) |>
  slice(1:2)
over_predicted
ames_train |>
  slice(over_predicted$.row) |>
  select(Gr_Liv_Area, Neighborhood, Year_Built, Bedroom_AbvGr, Full_Bath)


# Try on the other types of resamples
# v-fold with repeats
rf_res_r_folds <-
  rf_wflow |>
  fit_resamples(resamples = ames_r_folds, control = keep_pred)
rf_res_r_folds
rf_metrics_r_folds <- collect_metrics(rf_res_r_folds)

# bootstrapping
rf_res_boots <-
  rf_wflow |>
  fit_resamples(resamples = ames_boots, control = keep_pred)
rf_res_boots
rf_metrics_boots <- collect_metrics(rf_res_boots)

# validation set
val_res <- rf_wflow |> fit_resamples(resamples = val_set)
val_res
collect_metrics(val_res)

# parallel processing
parallel::detectCores(logical = FALSE) # physical cores
parallel::detectCores(logical = TRUE) # cores that can be used for independent processes


ls()
################################################################################

# Final code for this chapter from the book ----
library(tidymodels)
data(ames)
ames <- mutate(ames, Sale_Price = log10(Sale_Price))

set.seed(502)
ames_split <- initial_split(ames, prop = 0.80, strata = Sale_Price)
ames_train <- training(ames_split)
ames_test <- testing(ames_split)

ames_rec <-
  recipe(Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + Latitude + Longitude, data = ames_train) |>
  step_log(Gr_Liv_Area, base = 10) |>
  step_other(Neighborhood, threshold = 0.01) |>
  step_dummy(all_nominal_predictors()) |>
  step_interact(~ Gr_Liv_Area:starts_with("Bldg_Type_")) |>
  step_ns(Latitude, Longitude, deg_free = 20)

lm_model <- linear_reg() |> set_engine("lm")

lm_wflow <-
  workflow() |>
  add_model(lm_model) |>
  add_recipe(ames_rec)

lm_fit <- fit(lm_wflow, ames_train)

rf_model <-
  rand_forest(trees = 1000) |>
  set_engine("ranger") |>
  set_mode("regression")

rf_wflow <-
  workflow() |>
  add_formula(
    Sale_Price ~ Neighborhood + Gr_Liv_Area + Year_Built + Bldg_Type + Latitude + Longitude
  ) |>
  add_model(rf_model)

set.seed(1001)
ames_folds <- vfold_cv(ames_train, v = 10)

keep_pred <- control_resamples(save_pred = TRUE, save_workflow = TRUE)

set.seed(1003)
rf_res <- rf_wflow |> fit_resamples(resamples = ames_folds, control = keep_pred)
