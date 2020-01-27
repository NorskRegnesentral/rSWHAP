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


# @param n Number of samples to simulate
# @param mean Mean of untruncated normal distribution 
# @param sd Standard deviation of untruncated normal distribution
# @param a Lower truncation limit
# @param b Upper truncation limit
# @param approxn If TRUE approximately n samples will be generated. If FALSE exactly n samples will be generated at a higher computational cost.
rtnorm = function(n, mean, sd, a, b, approxn = TRUE) {
    q = c(a,b)
    p = pnorm(q, mean, sd)
    m = 5 # If we want this functionality, we must calcultate this so that the P(to few samples) is sufficiently small, say 1e-6.
    
    
    if(approxn) {
        nt = round(n/(p[2] - p[1]))
    } else {
        nt = 5*nt 
    }
    
    x = rnorm(nt, mean, sd)
    x = x[x > a & x <b]
    
    if(!approxn) {
        x = x[1:n]
    }
    return(x)
}


# @param n Number of samples to simulate
# @param mean Mean of untruncated normal distribution 
# @param sd Standard deviation of untruncated normal distribution
# @param lambda paramter in Box Cox transformation
# @param approxn If TRUE approximately n samples will be generated. If FALSE exactly n samples will be generated at a higher computational cost.
rpred = function(n, mean, sd, lambda, approxn = TRUE) {
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


# Example simulation from the model
n = 1000
mean = 2
sd = 2
lambda = 0.5
x = rpred(n, mean, sd, lambda, approxn = FALSE)
hist(x)

# Example computation of RMSE
mean = seq(5.1, 8, length.out = 10) # Mean of predictive distribution (in Box Cox world) in for ten temporal positions
sd = 1 # SD of predictive distribution (in Box Cox world) same in each temporal position
nobs = length(mean)
obs = rnorm(n = nobs, mean = 10, sd = 1) # Imagined observations orginal scale in each temporal position

# For every observation and prediction (temporal position) we must first estimate the expectation of the Box Cox distribution by simulation (lines 83 - 84). Next compute squared difference to observation. Som over position, finally take the quare root of the mean.
RMSE = 0
for(i in 1:nobs) {
    x = rpred(n, mean[i], sd, lambda, approxn = FALSE)
    Ex = mean(x)
    RMSE = RMSE + (obs[i] - Ex)^2
    RMSE = sqrt(RMSE/n)
}
