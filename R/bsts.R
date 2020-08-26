# Quick tutorial with bsts time series package
# http://www.unofficialgoogledatascience.com/2017/07/fitting-bayesian-structural-time-series.html

install.packages('bsts')

library(bsts)
data(iclaims)
str(initial.claims)
ss <- AddLocalLinearTrend(list(), initial.claims$iclaimsNSA)
ss <- AddSeasonal(ss, initial.claims$iclaimsNSA, nseasons=52)
model1 <- bsts(initial.claims$iclaimsNSA,
               state.specification=ss,
               niter=1000
               )
summary(model1)
plot(model1)
plot(model1, "components")
pred1 <- predict(model1, horizon=12)
plot(pred1, plot.original=156)

# Add a regression component
model2 <- bsts(iclaimsNSA ~ .,
               state.specification = ss,
               niter = 1000,
               data=initial.claims)

model3 <- bsts(iclaimsNSA ~ .,
               state.specification = ss,
               niter = 1000,
               data = initial.claims,
               expected.model.size=5  # extra argument for SpikeSlabPrior
               )

plot(model2, "components")
plot(model2, "coef")
pred2 <- predict(model2, horizon=12)
plot(pred2, plot.original=156)
bsts.prediction.errors(model2)

plot(model3, "components")
plot(model3, "coef")

# This plot compares residual errors ("mean absolute one step prediction errors")
# across multiple models
CompareBstsModels(list(
    "Model 1" = model1,
    "Model 2" = model2,
    "Model 3" = model3
    ),
    colors = c('black', 'green', 'blue')
)
# Adding regression on the Google search terms reduces the absolute error
# right around the spike in the 2008 recession.
# Models 2 & 3 do a similar job though.


