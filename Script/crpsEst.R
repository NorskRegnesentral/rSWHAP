# @param obst value to invert
# @param lambda parameter in Box Cox transformation
InvBoxCox = function(obst, lambda) {
    if(is.na(lambda)) {
        obs = NA
    } else {
        if(lambda == 0) {
            obs = exp(obst)
        } else {
            obs = (lambda*obst + 1)^(1/lambda)
        }
    }
    return(obs)
}


# @param n Approximately number of samples to generate
# @param mean Mean of untruncated normal distribution
# @param sd Standard deviation of untruncated normal distribution
# @param a Lower truncation limit
# @param b Upper truncation limit
rtnorm = function(n, mean, sd, a, b) {
    q = c(a,b)
    p = pnorm(q, mean, sd)
    nt = round(n/(p[2] - p[1]))

    x = rnorm(nt, mean, sd)
    x = x[x > a & x <b]
    return(x)
}


# @param n Number of samples to generate
# @param mean Mean of untruncated normal distribution
# @param sd Standard deviation of untruncated normal distribution
# @param lambda paramter in Box Cox transformation
rpred = function(n, mean, sd, lambda) {
    if(lambda < -1e-6) {
        a = -Inf
        b = -1/lambda
        xbc = rtnorm(n, mean, sd, a, b)
    } else if(lambda > 1e-6) {
        a = -1/lambda
        b = Inf
        xbc = rtnorm(n, mean, sd, a, b)
    } else {
        xbc = rnorm(n, mean, sd)
    }
    x = InvBoxCox(xbc, lambda)
    return(x)
}

rmseEst = function(obs,
                   mean,
                   sd,
                   lambda,
                   Nsamples = 10000,
                   print2screen = TRUE){

    nobs = length(obs)

    if(length(mean) != nobs | length(sd) != nobs | length(lambda) != nobs) {
        stop("Variables obs, mean, sd and lambda must contain the same number of elements.\n")
    }

    RMSE = 0
    for(i in 1:nobs) {
        x = rpred(Nsamples, mean[i], sd[i], lambda[i])
        Ex = mean(x)
        RMSE = RMSE + (obs[i] - Ex)^2
    }
    RMSE = sqrt(RMSE/nobs)
    return(RMSE)
}

crpsEst = function(obs, mean, sd, lambda, Nsamples = 10000) {
    nobs = length(obs)
    crps = rep(NA,nobs)

    #if(length(mean) != nobs | length(sd) != nobs | length(lambda) != nobs) {
    #    stop("Variables obs, mean, sd and lambda must contain the same number of elements.\n")
    #}

    for(i in 1:nobs) {
        EXz = 0
        EXXm = 0

        # Generate samples from Box Cox distribution
        x = rpred(Nsamples, mean[i], sd[i], lambda[i])
        EXz = mean(abs(x - obs))

        # Generate samples from Box Cox distribution
        x = rpred(Nsamples, mean[i], sd[i], lambda[i])
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


# Example computation of RMSE
# Syntetic data
mean = seq(5.1, 8, length.out = 10) # Mean of predictive distribution (in Box Cox world) in ten positions
obs = rnorm(n = length(mean), mean = 10, sd = 1) # Imagined observations on orginal scale in each position
nobs = length(obs)
lambda = rep(0.5, nobs)
sd = rep(1, nobs) # SD of predictive distribution (in Box Cox world) same in each position

# Estimating RMSE
rmseEst(obs, mean, sd, lambda, Nsamples = 1e5)

# Example error handling (sd with incorrect number of elements
rmseEst(obs, mean, sd[1], lambda, Nsamples = 1e5)

# Estimation of CRPS
crpsEst(obs, mean, sd, lambda, Nsamples = 1e5)

### TEST against C++ code ###
require(Rcpp)
sourceCpp("RMSE_CRPS_MEAN.cpp")
RMSE(obs, mean, sd[1], lambda[1], nsamples = 1e5)
mean(CRPS(obs, mean, sd[1], lambda[1], nsamples = 1e5))

