library(tidymodels)
tidymodels_prefer()

lm_model <- linear_reg() |>
    set_engine("lm")

lm_form_fit <- lm_model |>
    fit(Sale_Price ~ Longitude + Latitude, data = ames_train)

lm_xy_fit <- lm_model |>
    fit_xy(
        x = ames_train |> select(Longitude, Latitude),
        y = ames_train |> pull(Sale_Price)
    )

rf_model <- rand_forest(trees = 1000, min_n = 5) |>
    set_engine("ranger", verbose = TRUE) |>
    set_mode("regression")
rf_model |> translate()

# returns the fit object
lm_form_fit |> extract_fit_engine()

model_res <- lm_form_fit |>
    extract_fit_engine() |>
    summary()
model_res

param_est <- coef(model_res)
param_est
class(param_est)

# Get some tidy, consistently-formatted model results
broom::tidy(lm_form_fit)

ames_test_small <- ames_test |> slice(1:5)
predict(lm_form_fit, new_data = ames_test_small)

ames_test_small |>
    select(Sale_Price) |>
    # Add predictions
    bind_cols(predict(lm_form_fit, ames_test_small)) |>
    # Add 95% CIs
    bind_cols(predict(lm_form_fit, ames_test_small, type = "pred_int"))

tree_model <- decision_tree(min_n = 2) |>
    set_engine("rpart") |>
    set_mode("regression")

tree_fit <- tree_model |>
    fit(Sale_Price ~ Longitude + Latitude, data = ames_train)

ames_test_small |>
    select(Sale_Price) |>
    # Add predictions
    bind_cols(predict(tree_fit, ames_test_small))

# Get this nice template using the parsnip addin
parsnip::parsnip_addin()
decision_tree_C5.0_spec <-
    decision_tree(min_n = tune()) |>
    set_engine('C5.0')
