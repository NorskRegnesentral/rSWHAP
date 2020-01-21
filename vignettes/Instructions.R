## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3
)

## ----eval = FALSE-------------------------------------------------------------
#  devtools::install_github("NorskRegnesentral/rWHAP")

## ----eval=TRUE----------------------------------------------------------------
library(rWHAP)

## ----eval = TRUE--------------------------------------------------------------
data("ERAInterim")

## -----------------------------------------------------------------------------
summary(SWH)

## -----------------------------------------------------------------------------
summary(SLP)

## -----------------------------------------------------------------------------
summary(SLP.grad)

## -----------------------------------------------------------------------------
intercept.fourier = fourier(x = rep(1,length(time.all)))

## -----------------------------------------------------------------------------
training.test = split.data(years = years.all,
                           trainingPeriod = 2006:2014,
                           testPeriod = 2015)

## -----------------------------------------------------------------------------
pred.dist = getPreddistr(SWH = SWH,
          SLP = SLP, 
          SLP.grad = SLP.grad,
          latCell = 4, 
          longCell = 4,
          neig = 2,
          na.thresh = 500,
          latSWH = latitudeSWH,
          lonSWH = longitudeSWH,
          latSLP = latitudeSLP,
          longSLP = longitudeSLP,
          intercept.fourier = intercept.fourier,
          maxlag = 10)

# Extract the mean, standard deviation and the 
# estimated lambda parameter for the predictive
# distribution.
pred.mean = pred.dist$pred.mean
pred.sd = pred.dist$pred.sd
pred.lambda = pred.dist$pred.lambda

print(pred.dist$fits)

## -----------------------------------------------------------------------------
obs  <- SWH[4, 4, training.test[[2]]]

## -----------------------------------------------------------------------------
pit  <- pBoxCox(obs, pred.mean, pred.sd, pred.lambda)

## -----------------------------------------------------------------------------
mae = maeEst(obs = obs,
             mean = pred.mean,
             sd = pred.sd,
             lambda = pred.lambda)

## -----------------------------------------------------------------------------
plotPred(obs = obs,
         t.period = c(1:100),
         mean = pred.mean,
         sd = pred.sd,
         lambda = pred.lambda)

plotPred(obs = obs,
         t.period = c(1361:1460),
         mean = pred.mean,
         sd = pred.sd,
         lambda = pred.lambda)


## -----------------------------------------------------------------------------
rPlotPred(obs = obs,
         t.period = c(40:49),
         mean = pred.mean,
         sd = pred.sd,
         lambda = pred.lambda,
         n.random = 10)

