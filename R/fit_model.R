# Calculates the Fourier terms for modeling seasonality
# @param x Coefficients
# @param K Number of Fourier terms (default = 2)
# @param m Number of observations per period (default = 365.25*4)
# @return A matrix with 2xK Fourier terms
fourier = function(x, K = 2, m = 365.25*4) {
  n = length(x)
  idx = 1:n
  fourier.x = matrix(nrow = n, ncol = 2*K)
  coln = rep(NA, 2*K)
  for(k in 1:K) {
    fourier.x[,2*k-1] = sin(2*pi*k*idx/m)*x
    fourier.x[,2*k] = cos(2*pi*k*idx/m)*x
    coln[2*k-1] = paste("sin", k, sep = "")
    coln[2*k] = paste("cos", k, sep = "")
  }
  colnames(fourier.x) = coln
  colnames(fourier.x) = paste("intercept", colnames(fourier.x), sep = "_")
  return(fourier.x)
}


# Splitting the data in a traing set and a test set
# @param years A vector containg the year of each observation
# @param trainingPeriod A vector defining for which time period the model is trained on (default = 2006:2014)
# @param testPeriod Defining for which time period the model is tested on (default = 2015)
# @return A list containing the training set and the test set
split.data = function(years, trainingPeriod = 2006:2014, testPeriod = 2015) {

  training.test = list()
  training.test[[1]] = which(years %in% trainingPeriod)
  training.test[[2]] = which(years %in% testPeriod)

  return(training.test)
}

# Compute quantiles of the Box Cox distribution
# @param obs The probability
# @param mean The mean of the distribution
# @param sd The standard deviation of the distribution
# @param lambda The lambda used in the Box Cox transformation
# @return The quantiles of the Box Cox distribution
qBoxCox = function(p, mean, sd, lambda) {
  if(round(lambda,2) < 0) {
    q = (lambda*(sd*qnorm(p)+mean)+1)^(1/lambda)
  } else if(round(lambda,2) == 0) {
    q = exp(mean + sd*qnorm(p))
  } else { # lambda > 0
    T = (1/lambda + mean)/sd
    Vp = 1 - (1-p)*pnorm(T)
    q = (lambda*(sd*qnorm(Vp)+mean)+1)^(1/lambda)
  }
  return(q)
}

# Compute density of the Box Cox distribution
# @param obs The probability
# @param mean The mean of the distribution
# @param sd The standard deviation of the distribution
# @param lambda The lambda used in the Box Cox transformation
# @param log Logical value to switch between log or not
# @return The density values of the Box Cox distribution
dboxcox = function(obs, mean, sd, lambda, log = FALSE) {
  if(log == FALSE) {
    val = rep(0, length(obs))
    if(round(lambda,2) < 0) {
      val[obs>0] = obs[obs>0]^(lambda-1)*sd^(-1)*dnorm(((obs[obs>0]^lambda-1)/lambda-mean)/sd)*pnorm((-1/lambda-mean)/sd)^(-1)
    } else if(round(lambda,2) == 0) {
      val = dlnorm(obs, mean, sd, log = FALSE)
    } else {
      val[obs>0] = obs[obs>0]^(lambda-1)*sd^(-1)*dnorm(((obs[obs>0]^lambda-1)/lambda-mean)/sd)*pnorm((1/lambda+mean)/sd)^(-1)
    }
  } else {
    val = rep(-Inf, length(obs))
    if(round(lambda,2) < 0) {
      val[obs>0] = (lambda-1)*log(obs[obs>0]) - log(sd) + dnorm(((obs[obs>0]^lambda-1)/lambda-mean )/sd, log = TRUE) - log(pnorm( (-1/lambda - mean)/sd))
    } else if(round(lambda,2) == 0) {
      val = dlnorm(obs, mean, sd, log = TRUE)
    } else {
      val[obs>0] = (lambda-1)*log(obs[obs>0]) - log(sd) + dnorm(((obs[obs>0]^lambda-1)/lambda-mean)/sd, log = TRUE) - log(pnorm((1/lambda+mean)/sd))
    }
  }
  return(val)
}


# Compute the mean absolute error
# @param obs Observations
# @param mean The mean of the fitted distribution
# @param sd The standard deviation of the fitted distribution
# @param lambda The lambda used in the Box Cox transformation
# @return The mean absolute error of the Box Cox distribution
maeEst = function(obs,mean, sd, lambda) {

  median  <- qBoxCox(0.5, mean, sd, lambda)
  mae  <- mean(abs(obs - median))
  cat("\n The mean absolute error is:",mae,"\n")

  return(mae)
}

# Compute the log score
# @param obs Observations
# @param mean The mean of the fitted distribution
# @param sd The standard deviation of the fitted distribution
# @param lambda The lambda used in the Box Cox transformation
# @return The mean absolute error of the Box Cox distribution
logsEst = function(obs,mean, sd, lambda,log = TRUE) {

  logsD = dboxcox(obs = obs, mean = mean, sd = sd, lambda = lambda, log = log)
  logs = -mean(logsD)
  cat("\n The log score is:",logs,"\n")

  return(logs)
}

# Calculates the PIT values
# @param obs.bc Box Cox transformed observations
# @param mean The mean of the fitted distribution
# @param sd The standard deviation of the fitted distribution
# @param lambda The lambda used in the Box Cox transformation
# @return The PIT values normalized to the zero one intervall
ppredbc = function(obs, mean, sd, lambda) {
  if(round(lambda,2) < 0) {
    p = pnorm(obs = obs, mean = mean, sd = sd)/pnorm(-1/lambda, mean = mean, sd = sd)
  } else if(round(lambda,2) == 0) {
    p = pnorm(obs, mean = mean, sd = sd)
  } else { # lambda > 0
    p = pnorm(obs, mean = mean, sd = sd)/(1-pnorm(-1/lambda, mean = mean, sd = sd))
  }
  return(p)
}


# Computes the reliability index
# @param obs.bc Box Cox transformed observations
# @param mean The mean of the fitted distribution
# @param sd The standard deviation of the fitted distribution
# @param lambda The lambda used in the Box Cox transformation
# @param n.bins Number of bins used
# @return The reliability index (RIDX)
ribxEst = function(obs, mean, sd, lambda, n.bins = 20) {

  U = ppredbc(obs = obs,
             mean = mean,
             sd = sd,
             lambda = lambda)

  U = U[ !is.na(U) ]
  n = length(U)
  if(n > 0) {
    seps = seq(0, 1, length.out = n.bins+1)
    probs = rep(NA, n.bins)
    for(i in 1:n.bins) {
      probs[i] = sum( U >= seps[i] & U <= seps[i+1] )
    }
    probs = probs/n
    unif.prob = rep(1/n.bins, n.bins)
    RI = mean(abs(probs - unif.prob)*100)
  } else {
    RI = NA
  }

  cat("\n The reliability index is:",RI,"\n")

  return(RI)
}

# Computes the Continuous Ranked Probability Score (CRPS)
# @param obs.bc Box Cox transformed observations
# @param mean The mean of the fitted distribution
# @param sd The standard deviation of the fitted distribution
# @param lambda The lambda used in the Box Cox transformation
# @return The Continuous Ranked Probability Score (CRPS)
crpsEst = function(obs, mean, sd, lambda) {

  Rcpp::sourceCpp("RMSE_CRPS_MEAN.cpp")

  crpsEst = mean(CRPS(obs, mean, sd, lambda, n.bc.samp), na.rm = TRUE)
}

# Compute the root mean squared error
# @param obs Observations
# @param mean The mean of the fitted distribution
# @return The root mean squared error of the predictive distribution
rmseEst = function(obs,mean) {

  rmse  <- sqrt(mean(abs(obs - mean)))
  cat("\n The root mean squared error is:",mae,"\n")

  return(rmse)
}

# Calculates the continuous ranked probability score (CRPS)
# @param obs Observations
# @param mean The mean of the fitted distribution
# @return The root mean squared error of the predictive distribution
crpsEst = function(obs,mean) {

  crps  <- sqrt(mean(abs(obs - mean)))
  cat("\n The Continuous Ranked Probability Score (CRPS) is:",mae,"\n")

  return(crps)
}

# Performs Box-Cox transformation of the data. The transformation parameter (lambda) is chosen to minimize deviation from the normal distribution (minumum sum of squared skewness and squared kurtosis)
# @param obs The observations to be transformed
# @return The observations transformed
# @return The estimated lambda
BoxCoxLambda = function(obs) {
  lambda = seq(0.0, 1.0, 0.1)
  obs.trans = matrix(nrow = length(lambda), ncol = length(obs))
  normdev = rep(NA, length(lambda)) # Holds the amount of deviation from the normal distribution

  for(i in 1:length(lambda)) {
    if(lambda[i] == 0) {
      obs.trans[i,] = log(obs)
    } else {
      obs.trans[i,] = (obs^lambda[i]-1)/lambda[i]
    }
    normdev[i] = moments::skewness(obs.trans[i,],na.rm = TRUE)^2 + 0.25*(moments::kurtosis(obs.trans[i,],na.rm = TRUE))^2
  }
  return(list(data=obs.trans[which.min(normdev),],lambda = lambda[which.min(normdev)]))
}

# Plot the predictive distribution
# @param obs The observations to be transformed
# @param t.period The time period for plotting
# @param mean The mean value of the distribution
# @param sd The sd value of the distribution
# @param lambda The estimated lambda of the distribution
# @param graphDev Logical variable to create a OS specific graphic device
# @return A figure showing the predictive distribution
plotPred = function(obs,t.period, mean, sd,lambda, graphD=FALSE) {
  upper  <- qBoxCox(0.95, mean = mean, sd = sd, lambda = lambda)
  lower  <- qBoxCox(0.05, mean = mean, sd = sd, lambda = lambda)
  median  <- qBoxCox(0.5, mean, sd, lambda)

  t.upper  <- upper[t.period]
  t.lower  <- lower[t.period]


  if(graphD) graphDev(width = 7,height = 5)
  plot(t.period, obs[t.period], type="l",
       xlab="Time point in test period", ylab="SWH", ylim=c(0,15),
       main="")
  polygon(c(t.period, rev(t.period), t.period[1]),
          c(t.lower, rev(t.upper), t.lower[1]), col="gray90", border="NA")
  lines(t.period, obs[t.period])

  lines(t.period, median[t.period], col="gray50", lty=2)

  legend("topright", lty=c(1,2,1), lwd=c(1,1,4), col=c("black", "gray50", "gray90"),
         legend=c("Observation", "Predictive median", "90% prediction interval"))

}


# Plot random predictive trajectories for a selection of time points
# @param obs The observations to be transformed
# @param t.period The time period for plotting
# @param mean The mean value of the distribution
# @param sd The sd value of the distribution
# @param lambda The estimated lambda of the distribution
# @param n.random The number of random predictive trajectories
# @param graphDev Logical variable to create a OS specific graphic device
# @return A figure showing the random predictive trajectories
rPlotPred = function(obs,
                     t.period,
                     mean,
                     sd,
                     lambda,
                     n.random,
                     graphD = FALSE) {

  t.ind  = t.period
  n.period = length(t.ind)
  random.q  = array(NA, dim=c(n.random,n.period))
  for(i in 1:n.period) random.q[,i]  = qBoxCox(runif(n.random), mean[t.ind[i]], sd, lambda)

  if(graphD) graphDev(width = 7,height = 5)
  plot(t.ind, random.q[1,], type="l", col="gray50",
       xlab="Time point in test period", ylab="SWH", ylim=c(5,12),
       main="Random predictive trajectories")
  for(i in 2:n.random) lines(t.ind, random.q[i,], col="gray50")
  lines(t.ind, obs[t.ind], col="black", lwd=2)

}

# Plot random correlation for a selection of time points
# @param obs The observations to be transformed
# @param t.period The time period for plotting
# @param mean The mean value of the distribution
# @param sd The sd value of the distribution
# @param lambda The estimated lambda of the distribution
# @param n.random The number of random predictive trajectories
# @param training.test The training and the test period
# @param SWHobs SWH observations for a particular location in the test period
# @param graphDev Logical variable to create a OS specific graphic device
# @return A figure showing the random correlations
rCorr = function(obs,
                 t.period,
                 mean,
                 sd,
                 lambda,
                 n.random,
                 training.test,
                 SWHobs,
                 graphD = FALSE) {


  n.period = length(t.period)
  t.ind  = t.period

  random.q  = array(NA, dim=c(n.random,n.period))
  for(i in 1:n.period) random.q[,i]  = qBoxCox(runif(n.random), pred.mean[t.ind[i]], pred.sd, pred.lambda)

  sample.q  <- array(NA, dim=c(n.random,n.period))
  for(i in 1:n.period) sample.q[,i]  <- rank(random.q[,i])
  sample.q

  nTest  <- length(training.test[[2]])
  h.ind  <- training.test[[2]][(nTest-(n.random*n.period-1)):nTest]
  hist.obs  <- t(array(SWHobs[h.ind], dim=c(n.period,n.random)))
  hist.q  <- array(NA, dim=c(n.random,n.period))
  for(i in 1:n.period) hist.q[,i]  = rank(hist.obs[,i])
  hist.q
  sort.q  <- random.q
  for(i in 1:n.period) sort.q[,i] = sort(random.q[,i])[hist.q[,i]]


  if (graphD) graphDev(width = 7,height = 5)
  plot(t.ind, sort.q[1,], type="l", col="gray50",
       xlab="Time point in test period", ylab="SWH", ylim=c(5,12),main="")
  for(i in 2:n.random) lines(t.ind, sort.q[i,], col="gray50")
  lines(t.ind, obs[t.ind], col="black", lwd=2)
}

# Prephare graphical device for a given OS
# @param width Width of the graphical window
# @param height Height of the grapchical window
# @return The an OS dependent graphical window
graphDev = function(width = 7,height = 5) {

  system = Sys.info()
  if(system['sysname']=="Windows"){
    windows(width = 7,height = 5)
  }

  if(system['sysname']=="Linux"){
    X11(width = 7,height = 5)
  }

  if(system['sysname']=="Darwin"){
    quartz("",width = 7,height = 5)
  }
}

# Compute PIT values of the Box Cox distribution
# @param obs The observations to be transformed
# @param mean The mean of the predictive distribution
# @param sd The standard deviation of the predictive distribution
# @param lambda The estimated Box Cox lambda parameter used in the transformation of the predictive distribution
# @param plothist Set to TRUE if plotting the histogram of the PIT values
# @param graphDev Logical variable to create a OS specific graphic device
# @return The PIT values
pBoxCox = function(x, mean, sd, lambda,plothist = TRUE,graphD = FALSE) {
  g.x  <- BoxCoxLambdaKnown(x, lambda)
  if(round(lambda,2) < 0) {
    p = pnorm((g.x-mean)/sd) / pnorm((-1/lambda - mean)/sd)
  } else if(round(lambda,2) == 0) {
    p = pnorm((g.x-mean)/sd)
  } else { # lambda > 0
    p = (pnorm((g.x-mean)/sd) - pnorm((-1/lambda - mean)/sd)) / pnorm((1/lambda + mean)/sd)
  }

  if(graphD) graphDev(width = 7,height = 5)
  hist(p, freq=FALSE, nclass=10, col="gray", xlab="PIT value",main = "PIT histogram")
  abline(a=1, b=0, lty=2,col = "red")
  return(p)
}

# Performs Box-Cox transformation with a known lambda parameter
# @param obs The observations to be transformed
# @param lambda The lambda parameter to be used in the transformation
# @return The observations transformed
BoxCoxLambdaKnown = function(obs, lambda) {
  if(round(lambda,3) == 0) {
    obs = log(obs)
  } else {
    obs = (obs^lambda - 1)/lambda
  }
  return(obs)
}

# Calculates the AR lags
# @param n.training Size of training data
# @param n.test Size of training data
# @param maxlag Order of AR process
# @return AR terms for various lags used in model fitting
makeARterms = function(SWH.bc.standard.training,
                       SWH.bc.fourier.training,
                       SWH.bc.standard.test,
                       SWH.bc.fourier.test,
                       n.training = 100,
                       n.test = 10,
                       maxlag = 10) {

  SWH.bc.ar.training.names = rep("",maxlag)
  SWH.bc.fourier.ar.training.names = rep("",maxlag)
  SWH.bc.ar.test.names = rep("",maxlag)
  SWH.bc.fourier.ar.test.names = rep("",maxlag)

  SWH.bc.ar.training.list = list()
  SWH.bc.fourier.ar.training.list = list()
  SWH.bc.ar.test.list = list()
  SWH.bc.fourier.ar.test.list = list()

  for(lag in 1:maxlag){
    #assign(paste("SWH.bc.ar.training.m",lag,sep = ""),c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)]))
    SWH.bc.ar.training.names[lag] = paste("SWH.bc.ar.training.m",lag,sep = "")
    SWH.bc.ar.training.list[[lag]] = c(rep(SWH.bc.standard.training[1],lag), SWH.bc.standard.training[1:(n.training-lag)])

    #assign(paste("SWH.bc.fourier.ar.training.m",lag,sep=""),rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),]))
    SWH.bc.fourier.ar.training.names[lag] = paste("SWH.bc.fourier.ar.training.m",lag,sep = "")
    SWH.bc.fourier.ar.training.list[[lag]] = rbind(matrix(rep(SWH.bc.fourier.training[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.training[1:(n.training-lag),])

    #assign(paste("SWH.bc.ar.test.m",lag,sep = ""), c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)]))
    SWH.bc.ar.test.names[lag] = paste("SWH.bc.ar.test.m",lag,sep = "")
    SWH.bc.ar.test.list[[lag]] = c(rep(SWH.bc.standard.test[1],lag), SWH.bc.standard.test[1:(n.test-lag)])

    #assign(paste("SWH.bc.fourier.ar.test.m",lag,sep = ""), rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),]))
    SWH.bc.fourier.ar.test.names[lag] = paste("SWH.bc.fourier.ar.test.m",lag,sep = "")
    SWH.bc.fourier.ar.test.list[[lag]] = rbind(matrix(rep(SWH.bc.fourier.test[1,],lag), nrow = lag, byrow = TRUE), SWH.bc.fourier.test[1:(n.test-lag),])

  }

  names(SWH.bc.ar.training.list) = SWH.bc.ar.training.names
  names(SWH.bc.fourier.ar.training.list) = SWH.bc.fourier.ar.training.names
  names(SWH.bc.ar.test.list) = SWH.bc.ar.test.names
  names(SWH.bc.fourier.ar.test.list) = SWH.bc.fourier.ar.test.names



  ARterms = list()
  ARterms$SWH.bc.ar.training.list = SWH.bc.ar.training.list
  ARterms$SWH.bc.fourier.ar.training.list = SWH.bc.fourier.ar.training.list
  ARterms$SWH.bc.ar.test.list = SWH.bc.ar.test.list
  ARterms$SWH.bc.fourier.ar.test.list = SWH.bc.fourier.ar.test.list

  return(ARterms)
}



# Training the model and generating the predictive distributions
# @param SWH The SWH data
# @param SLP The SLP data
# @param SLP.grad The SLP gradient data
# @param latCell Define the latitude cell of interest (default = 4)
# @param longCell Define the longitude cell of interest (default = 4)
# @param training.test The training and the test data
# @param neig Define the number of spatial neighbours (default = 2)
# @param na.thresh Define how many missing values is "good enough" (default = 500)
# @param latSWH A vector containing the SWH latitudes
# @param lonSWH A vector containing the SWH longitudes
# @param latSLP A vector containing the SLP latitudes
# @param lonSLP A vector containing the SLO longitudes
# @param intercept.fourier Seasonal Fourier coefficients
# @param maxlag Define the maximum order of AR processes to be constructed
# @return A list of predictive means, predictive spread and lambda for Box-Cox transformation
getPreddistr = function(SWH = NA,
                     SLP = NA,
                     SLP.grad = NA,
                     latCell = 4,
                     longCell = 4,
                     neig = 2,
                     na.thresh = 500,
                     latSWH = NA,
                     lonSWH = NA,
                     latSLP = NA,
                     longSLP = NA,
                     intercept.fourier = NA,
                     maxlag = 10) {

  pred.mean = rep(NA, length(training.test[[2]]))
  pred.sd = NA
  pred.lambda = NA

  ## Divide in training and test set
  idx.training = training.test[[1]]
  idx.test = training.test[[2]]

  idx.rank = 0

  ## Make sure covariates are from the correct location
  idx.longSLP = which(longitudeSLP == longitudeSWH[longCell])
  idx.latSLP = which(latitudeSLP == latitudeSWH[latCell])

  ## Only continue with analysis if sufficient available data
  if(sum(is.na(SWH[longCell, latCell,])) < na.thresh &
     sum(is.na( SLP[idx.longSLP, idx.latSLP, ])) < na.thresh &
     sum(is.na( SLP.grad[idx.longSLP, idx.latSLP, ])) < na.thresh) {

    ## Remove missing data
    not.NA = which(!is.na(SWH[longCell, latCell,]) &
                     !is.na( SLP[idx.longSLP, idx.latSLP, ]) &
                     !is.na( SLP.grad[idx.longSLP, idx.latSLP, ]))
    idx.training = idx.training[idx.training %in% not.NA]
    n.training = length(idx.training)
    idx.test = idx.test[idx.test %in% not.NA]
    n.test = length(idx.test)

    ## Estimate lambda and Box-Cox transform SWH in training period
    SWH.bc.dummy.training = BoxCoxLambda(SWH[longCell, latCell,idx.training])
    SWH.bc.lambda.training = SWH.bc.dummy.training$lambda
    SWH.bc.training = SWH.bc.dummy.training$data

    ## Box-Cox transform SWH in test period with same lambda
    SWH.bc.test = BoxCoxLambdaKnown(SWH[longCell, latCell,idx.test], SWH.bc.lambda.training)

    ## Estimate lambda and Box-Cox transform SLP.grad in training period
    SLP.grad.bc.dummy.training = BoxCoxLambda(SLP.grad[idx.longSLP, idx.latSLP,idx.training])
    SLP.grad.bc.lambda.training = SLP.grad.bc.dummy.training$lambda
    SLP.grad.bc.training = SLP.grad.bc.dummy.training$data
    ## Box-Cox transform SLP.grad in test period with same lambda
    SLP.grad.bc.test = BoxCoxLambdaKnown(SLP.grad[idx.longSLP, idx.latSLP,idx.test], SLP.grad.bc.lambda.training)

    ## Standardizing the variables
    SWH.bc.mean.training = mean(SWH.bc.training)
    SWH.bc.sd.training = sd(SWH.bc.training)
    SWH.bc.standard.training = (SWH.bc.training - SWH.bc.mean.training)/
      SWH.bc.sd.training
    SWH.bc.standard.test = (SWH.bc.test - SWH.bc.mean.training)/
      SWH.bc.sd.training

    SLP.mean.training = mean( SLP[idx.longSLP, idx.latSLP, idx.training] )
    SLP.sd.training = sd( SLP[idx.longSLP, idx.latSLP, idx.training] )
    SLP.standard.training = (SLP[idx.longSLP, idx.latSLP, idx.training] - SLP.mean.training)/
      SLP.sd.training
    SLP.standard.test = (SLP[idx.longSLP, idx.latSLP, idx.test] - SLP.mean.training)/
      SLP.sd.training

    SLP.grad.bc.mean.training = mean(SLP.grad.bc.training)
    SLP.grad.bc.sd.training = sd(SLP.grad.bc.training)
    SLP.grad.bc.standard.training = (SLP.grad.bc.training - SLP.grad.bc.mean.training)/
      SLP.grad.bc.sd.training
    SLP.grad.bc.standard.test = (SLP.grad.bc.test - SLP.grad.bc.mean.training)/
      SLP.grad.bc.sd.training

    ## Compute Fourier periodics of the explanatory variables
    fourier.training = intercept.fourier[idx.training,]
    fourier.test = intercept.fourier[idx.test,]

    SLP.fourier.training = fourier.training*SLP.standard.training
    SLP.fourier.test = fourier.test*SLP.standard.test

    SLP.grad.bc.fourier.training = fourier.training*SLP.grad.bc.standard.training
    SLP.grad.bc.fourier.test = fourier.test*SLP.grad.bc.standard.test

    SWH.bc.fourier.training = fourier.training*SWH.bc.standard.training
    SWH.bc.fourier.test = fourier.test*SWH.bc.standard.test

    # Create AR terms for model selection
    ARterms = makeARterms(SWH.bc.standard.training = SWH.bc.standard.training,
                          SWH.bc.fourier.training = SWH.bc.fourier.training,
                          SWH.bc.standard.test = SWH.bc.standard.test,
                          SWH.bc.fourier.test = SWH.bc.fourier.test,
                          n.training = n.training,
                          n.test = n.test,
                          maxlag = maxlag)

    SWH.bc.ar.training.list = ARterms$SWH.bc.ar.training.list
    SWH.bc.fourier.ar.training.list = ARterms$SWH.bc.fourier.ar.training.list
    SWH.bc.ar.test.list = ARterms$SWH.bc.ar.test.list
    SWH.bc.fourier.ar.test.list = ARterms$SWH.bc.fourier.ar.test.list

    #Extract variables from the lists
    list2env(setNames(SWH.bc.ar.training.list,paste0("SWH.bc.ar.training.m",seq_along(SWH.bc.ar.training.list))), envir = parent.frame())
    list2env(setNames(SWH.bc.fourier.ar.training.list,paste0("SWH.bc.fourier.ar.training.m",seq_along(SWH.bc.fourier.ar.training.list))), envir = parent.frame())
    list2env(setNames(SWH.bc.ar.test.list,paste0("SWH.bc.ar.test.m",seq_along(SWH.bc.ar.test.list))), envir = parent.frame())
    list2env(setNames(SWH.bc.fourier.ar.test.list,paste0("SWH.bc.fourier.ar.test.m",seq_along(SWH.bc.fourier.ar.test.list))), envir = parent.frame())

        ## Vanem&Walker spatial model LASSO #####

    ## Compute mean, max og min of the neighborhood of the current point
    SLP.spatmax = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                             max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = TRUE)
    SLP.spatmin = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                             max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = T)
    SLP.spatmean = apply(SLP[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                              max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = TRUE)

    SLP.grad.spatmax = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                       max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, max, na.rm = TRUE )
    SLP.grad.spatmin = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                       max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, min, na.rm = TRUE)
    SLP.grad.spatmean = apply(SLP.grad[ max(idx.longSLP-neig,1):min(idx.longSLP+neig, dim(SLP)[1]),
                                        max(idx.latSLP-neig,1):min(idx.latSLP+neig, dim(SLP)[2]),], 3, mean, na.rm = TRUE)

    ## Split in test and training
    SLP.spatmax.training = SLP.spatmax[idx.training]
    SLP.spatmin.training = SLP.spatmin[idx.training]
    SLP.spatmean.training = SLP.spatmean[idx.training]

    SLP.grad.spatmax.training = SLP.grad.spatmax[idx.training]
    SLP.grad.spatmin.training = SLP.grad.spatmin[idx.training]
    SLP.grad.spatmean.training = SLP.grad.spatmean[idx.training]

    SLP.spatmax.test = SLP.spatmax[idx.test]
    SLP.spatmin.test = SLP.spatmin[idx.test]
    SLP.spatmean.test = SLP.spatmean[idx.test]

    SLP.grad.spatmax.test = SLP.grad.spatmax[idx.test]
    SLP.grad.spatmin.test = SLP.grad.spatmin[idx.test]
    SLP.grad.spatmean.test = SLP.grad.spatmean[idx.test]

    ## Standardize SLP spatial variables
    SLP.spatmax.mean.training = mean(SLP.spatmax.training)
    SLP.spatmax.sd.training = sd(SLP.spatmax.training)
    SLP.spatmax.standard.training = (SLP.spatmax.training - SLP.spatmax.mean.training)/
      SLP.spatmax.sd.training
    SLP.spatmax.standard.test = (SLP.spatmax.test - SLP.spatmax.mean.training)/
      SLP.spatmax.sd.training

    SLP.spatmin.mean.training = mean(SLP.spatmin.training)
    SLP.spatmin.sd.training = sd(SLP.spatmin.training)
    SLP.spatmin.standard.training = (SLP.spatmin.training - SLP.spatmin.mean.training)/
      SLP.spatmin.sd.training
    SLP.spatmin.standard.test = (SLP.spatmin.test - SLP.spatmin.mean.training)/
      SLP.spatmin.sd.training

    SLP.spatmean.mean.training = mean(SLP.spatmean.training)
    SLP.spatmean.sd.training = sd(SLP.spatmean.training)
    SLP.spatmean.standard.training = (SLP.spatmean.training - SLP.spatmean.mean.training)/
      SLP.spatmean.sd.training
    SLP.spatmean.standard.test = (SLP.spatmean.test - SLP.spatmean.mean.training)/
      SLP.spatmean.sd.training

    ## Box-Cox transform and standardize SLP.grad
    SLP.grad.spatmax.bc.dummy = BoxCoxLambda(SLP.grad.spatmax.training)
    SLP.grad.spatmax.bc.lambda.training = SLP.grad.spatmax.bc.dummy$lambda
    SLP.grad.spatmax.bc.training = SLP.grad.spatmax.bc.dummy$data
    SLP.grad.spatmax.bc.test = BoxCoxLambdaKnown(SLP.grad.spatmax.test, SLP.grad.spatmax.bc.lambda.training)

    SLP.grad.spatmax.bc.mean.training = mean(SLP.grad.spatmax.bc.training)
    SLP.grad.spatmax.bc.sd.training = sd(SLP.grad.spatmax.bc.training)
    SLP.grad.spatmax.bc.standard.training = (SLP.grad.spatmax.bc.training - SLP.grad.spatmax.bc.mean.training)/
      SLP.grad.spatmax.bc.sd.training
    SLP.grad.spatmax.bc.standard.test = (SLP.grad.spatmax.bc.test - SLP.grad.spatmax.bc.mean.training)/
      SLP.grad.spatmax.bc.sd.training

    SLP.grad.spatmin.bc.dummy = BoxCoxLambda(SLP.grad.spatmin.training)
    SLP.grad.spatmin.bc.lambda.training = SLP.grad.spatmin.bc.dummy$lambda
    SLP.grad.spatmin.bc.training = SLP.grad.spatmin.bc.dummy$data
    SLP.grad.spatmin.bc.test = BoxCoxLambdaKnown(SLP.grad.spatmin.test, SLP.grad.spatmin.bc.lambda.training)

    SLP.grad.spatmin.bc.mean.training = mean(SLP.grad.spatmin.bc.training)
    SLP.grad.spatmin.bc.sd.training = sd(SLP.grad.spatmin.bc.training)
    SLP.grad.spatmin.bc.standard.training = (SLP.grad.spatmin.bc.training - SLP.grad.spatmin.bc.mean.training)/
      SLP.grad.spatmin.bc.sd.training
    SLP.grad.spatmin.bc.standard.test = (SLP.grad.spatmin.bc.test - SLP.grad.spatmin.bc.mean.training)/
      SLP.grad.spatmin.bc.sd.training

    SLP.grad.spatmean.bc.dummy = BoxCoxLambda(SLP.grad.spatmean.training)
    SLP.grad.spatmean.bc.lambda.training = SLP.grad.spatmean.bc.dummy$lambda
    SLP.grad.spatmean.bc.training = SLP.grad.spatmean.bc.dummy$data
    SLP.grad.spatmean.bc.test = BoxCoxLambdaKnown(SLP.grad.spatmean.test, SLP.grad.spatmean.bc.lambda.training)

    SLP.grad.spatmean.bc.mean.training = mean(SLP.grad.spatmean.bc.training)
    SLP.grad.spatmean.bc.sd.training = sd(SLP.grad.spatmean.bc.training)
    SLP.grad.spatmean.bc.standard.training = (SLP.grad.spatmean.bc.training - SLP.grad.spatmean.bc.mean.training)/
      SLP.grad.spatmean.bc.sd.training
    SLP.grad.spatmean.bc.standard.test = (SLP.grad.spatmean.bc.test - SLP.grad.spatmean.bc.mean.training)/
      SLP.grad.spatmean.bc.sd.training

    ## Compute Fourier predictors for seasonally varying coefficients
    SLP.spatmax.fourier.training = fourier.training*SLP.spatmax.standard.training
    SLP.spatmax.fourier.test = fourier.test*SLP.spatmax.standard.test

    SLP.spatmin.fourier.training = fourier.training*SLP.spatmin.standard.training
    SLP.spatmin.fourier.test = fourier.test*SLP.spatmin.standard.test

    SLP.spatmean.fourier.training = fourier.training*SLP.spatmean.standard.training
    SLP.spatmean.fourier.test = fourier.test*SLP.spatmean.standard.test

    SLP.grad.spatmax.bc.fourier.training = fourier.training*SLP.grad.spatmax.bc.standard.training
    SLP.grad.spatmax.bc.fourier.test = fourier.test*SLP.grad.spatmax.bc.standard.test

    SLP.grad.spatmin.bc.fourier.training = fourier.training*SLP.grad.spatmin.bc.standard.training
    SLP.grad.spatmin.bc.fourier.test = fourier.test*SLP.grad.spatmin.bc.standard.test

    SLP.grad.spatmean.bc.fourier.training = fourier.training*SLP.grad.spatmean.bc.standard.training
    SLP.grad.spatmean.bc.fourier.test = fourier.test*SLP.grad.spatmean.bc.standard.test

    ## Create full set of spatial predictors for the training and test sets
    predictors.training.spatial = as.data.frame(cbind(SLP.spatmax.standard.training,
                                                      SLP.spatmin.standard.training,
                                                      SLP.spatmean.standard.training,
                                                      SLP.grad.spatmax.bc.standard.training,
                                                      SLP.grad.spatmin.bc.standard.training,
                                                      SLP.grad.spatmean.bc.standard.training,
                                                      SLP.spatmax.fourier.training,
                                                      SLP.spatmin.fourier.training,
                                                      SLP.spatmean.fourier.training,
                                                      SLP.grad.spatmax.bc.fourier.training,
                                                      SLP.grad.spatmin.bc.fourier.training,
                                                      SLP.grad.spatmean.bc.fourier.training))

    predictors.test.spatial = as.data.frame(cbind(SLP.spatmax.standard.test,
                                                  SLP.spatmin.standard.test,
                                                  SLP.spatmean.standard.test,
                                                  SLP.grad.spatmax.bc.standard.test,
                                                  SLP.grad.spatmin.bc.standard.test,
                                                  SLP.grad.spatmean.bc.standard.test,
                                                  SLP.spatmax.fourier.test,
                                                  SLP.spatmin.fourier.test,
                                                  SLP.spatmean.fourier.test,
                                                  SLP.grad.spatmax.bc.fourier.test,
                                                  SLP.grad.spatmin.bc.fourier.test,
                                                  SLP.grad.spatmean.bc.fourier.test))

    ## Create set of local predictors with AR order 5
    predictors.training = as.data.frame( cbind(predictors.training.spatial,
                                               intercept.fourier[idx.training,],
                                               SLP.standard.training,
                                               SLP.grad.bc.standard.training,
                                               SLP.fourier.training,
                                               SLP.grad.bc.fourier.training,
                                               SWH.bc.ar.training.m1,
                                               SWH.bc.fourier.ar.training.m1,
                                               SWH.bc.ar.training.m2,
                                               SWH.bc.fourier.ar.training.m2,
                                               SWH.bc.ar.training.m3,
                                               SWH.bc.fourier.ar.training.m3,
                                               SWH.bc.ar.training.m4,
                                               SWH.bc.fourier.ar.training.m4,
                                               SWH.bc.ar.training.m5,
                                               SWH.bc.fourier.ar.training.m5) )

    predictors.test = as.data.frame( cbind(predictors.test.spatial,
                                           intercept.fourier[idx.test,],
                                           SLP.standard.test,
                                           SLP.grad.bc.standard.test,
                                           SLP.fourier.test,
                                           SLP.grad.bc.fourier.test,
                                           SWH.bc.ar.test.m1,
                                           SWH.bc.fourier.ar.test.m1,
                                           SWH.bc.ar.test.m2,
                                           SWH.bc.fourier.ar.test.m2,
                                           SWH.bc.ar.test.m3,
                                           SWH.bc.fourier.ar.test.m3,
                                           SWH.bc.ar.test.m4,
                                           SWH.bc.fourier.ar.test.m4,
                                           SWH.bc.ar.test.m5,
                                           SWH.bc.fourier.ar.test.m5) )

    colnames(predictors.training) = paste("V", 1:dim(predictors.training)[2], sep = "")
    colnames(predictors.test) = paste("V", 1:dim(predictors.test)[2], sep = "")
    cat("Total number of potential predictors used in LASSO selection:",dim(predictors.training)[2], "\n")

    ## Start with LASSO selection
    cv = glmnet::cv.glmnet(as.matrix(predictors.training), SWH.bc.standard.training, family = "gaussian", alpha = 1, nfold = 10)
    lasso = glmnet::glmnet(as.matrix(predictors.training), SWH.bc.standard.training, alpha = 1)
    minindex = which.min(abs(lasso$lambda - cv$lambda.min))
    beta = lasso$beta[,minindex]
    predictors.training = predictors.training[, which( abs(beta) > 1e-6 )]
    predictors.test = predictors.test[, which( abs(beta) > 1e-6 )]
    cat("Number of predictors selected by LASSO:",dim(predictors.training)[2], "\n")

    cat("Fitting linear model with the variables from the LASSO selection...\n")
    fit <- lm(SWH.bc.standard.training ~ ., data = predictors.training)
    fits = summary(fit)

    ## predict
    SWH.bc.standard.pred <- predict(object = fit, newdata = predictors.test) #lm
    SWH.bc.standard.pred.se <- fits$sigma

    ## Descale and detransform before computing different raknkings
    SWH.bc.pred = SWH.bc.standard.pred*SWH.bc.sd.training + SWH.bc.mean.training
    SWH.bc.pred.se = SWH.bc.standard.pred.se*SWH.bc.sd.training
    SWH.pred = forecast::InvBoxCox(SWH.bc.pred, SWH.bc.lambda.training) #detransform

    pred.mean[idx.test - length(training.test[[1]])] = SWH.bc.pred
    pred.sd = SWH.bc.pred.se
    pred.lambda = SWH.bc.lambda.training
  }
  return(list(pred.SWH=SWH.pred, pred.mean=pred.mean, pred.sd=pred.sd, pred.lambda=pred.lambda, fits=fits))

}


