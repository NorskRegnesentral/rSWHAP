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

x = rpred(n = 1000, mean = 2, sd = 2, lambda = 0.5, approxn = FALSE)
x
hist(x)
