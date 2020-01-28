crpsEst = function(obs, mean, sd, lambda,Nsamples = 10000) {

  nobs = length(mean)
  crps = rep(NA,nobs)

  for(i in 1:nobs) {
    EXz = 0
    EXXm = 0

    # Generate samples from Box Cox distribution
    if(lambda < -1e-6) {
      a = -infty
      b = -1/lambda
      xbc = rtnorm(nsamples, mean[i], sd, a, b)
    } else if(lambda > 1e-6) {
      a = -1/lambda
      b = infty
      xbc = rtnorm(Nsamples, mean[i], sd, a, b)
    } else {
      xbc = rnorm(Nsamples, mean[i], sd)
    }
    x = InvBoxCox(xbc, lambda)

    for(j in 1:Nsamples) {
      EXz = EXz + abs(x[j] - obs[i])
    }
    EXz = EXz/Nsamples

    # Generate samples from Box Cox distribution
    if(lambda < -1e-6) {
      a = -infty
      b = -1/lambda
      xbc = rtnorm(Nsamples, mean[i], sd, a, b)
    } else if(lambda > 1e-6) {
      a = -1/lambda
      b = infty
      xbc = rtnorm(Nsamples, mean[i], sd, a, b)
    } else {
      xbc = rnorm(Nsamples, mean[i], sd)
    }
    x = InvBoxCox(xbc, lambda)

    nsamp = 0.0;
    for(j in 0:(Nsamples/2-1)) {
      EXXm = EXXm + abs(x[2*j+1] - x[2*j+2])
      nsamp = nsamp + 1.0
    }
    EXXm = EXXm/nsamp

    crps[i] = EXz - 0.5*EXXm;
  }

  crps = mean(crps)

  return(crps);

}
