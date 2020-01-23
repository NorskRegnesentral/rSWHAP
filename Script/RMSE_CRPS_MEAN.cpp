#include <Rcpp.h> 
using namespace Rcpp;

NumericVector InvBoxCox(NumericVector xbc, double lambda) {
    int i;    
    NumericVector x(xbc.size());
    for(i = 0; i<xbc.size(); i++) {
        if(lambda < -1e-6 || lambda > 1e-6) {
            x[i] = pow(lambda*xbc[i]+1,1/lambda);
        } else {
            x[i] = exp(xbc[i]);
        }
    }
    return(x);
}

// Generate approximately N samples from rtnorm
NumericVector rtnorm(int N, double mean, double sd, double a, double b) {
    // Compute number of samples needed from rnorm to get approximately N samples from tnorm
    NumericVector q(2);
    q[0] = a;
    q[1] = b;    
    NumericVector p = pnorm(q, mean, sd);
    int Nt = nearbyint(N/(p[1] - p[0])); 
    
    NumericVector samples = rnorm(Nt, mean, sd);
    return( samples[samples>a & samples<b] );
}

// [[Rcpp::export]]
NumericVector meanboxcox(NumericVector mean, double sd, double lambda, int nsamples) {    
    int i;
    int n = mean.size(); // mean - prediction mean on Box Cox scale
    double a, b, rmse = 0.0, infty = 1e308; 

    NumericVector xbc; // Samples from the Box Cox world (truncated normal)
    NumericVector x; // Samples from the org world
    NumericVector Ep(n); // Ep - Estimation of the expectation of the Box Cox distribution
    for(i=0; i<n; i++) {
        if(mean[i] == mean[i]) { // mean[i] is not NAN
            // Generate samples from Box Cox distribution
            if(lambda < -1e-6) {
                a = -infty;
                b = -1/lambda;
                xbc = rtnorm(nsamples, mean[i], sd, a, b);
            } else if(lambda > 1e-6) {
                a = -1/lambda;
                b = infty;
                xbc = rtnorm(nsamples, mean[i], sd, a, b);            
            } else {
                xbc = rnorm(nsamples, mean[i], sd);
            }
            x = InvBoxCox(xbc, lambda);
            Ep[i] = sum(x)/x.size(); // For some reason mean(x) don't work..
        } else {
            Ep[i] = NAN;
        }
    }
    return(Ep);
}    

// [[Rcpp::export]]
double RMSE(NumericVector w, NumericVector mean, double sd, double lambda, int nsamples) {
    int i;
    int n = w.size(); // w - the observations
    double Ep, a, b, rmse = 0.0, infty = 1e308; // Ep - Estimation of the expectation of the Box Cox distribution
    //NumericVector rmse(n);
    NumericVector xbc; // Samples from the Box Cox world (truncated normal)
    NumericVector x; // Samples from the org world
    for(i=0; i<n; i++) {
        // Generate samples from Box Cox distribution
        if(lambda < -1e-6) {
            a = -infty;
            b = -1/lambda;
            xbc = rtnorm(nsamples, mean[i], sd, a, b);
        } else if(lambda > 1e-6) {
            a = -1/lambda;
            b = infty;
            xbc = rtnorm(nsamples, mean[i], sd, a, b);            
        } else {
            xbc = rnorm(nsamples, mean[i], sd);
        }
        x = InvBoxCox(xbc, lambda);
        Ep = sum(x)/x.size(); // For some reason mean(x) don't work..
        //Rprintf("Ep = %f, mean[i] = %f, w[i] = %f \n", Ep, mean[i], w[i]);

        rmse += pow(w[i] - Ep,2);
    }
    return(sqrt(rmse/n));
}



// [[Rcpp::export]]
NumericVector CRPS(NumericVector z, NumericVector mean, double sd, double lambda, int nsamples) {
    int i, j;
    int n = z.size(); // z - the observations
    double EXz, EXXm, a, b, infty = 1e308, nsamp; // Ep - Estimation of the expectation of the Box Cox distribution
    NumericVector crps(n);
    NumericVector xbc; // Samples from the Box Cox world (truncated normal)
    NumericVector x; // Samples from the org world
    for(i=0; i<n; i++) {
        EXz = 0;
        EXXm = 0;
        
        // Generate samples from Box Cox distribution
        if(lambda < -1e-6) {
            a = -infty;
            b = -1/lambda;
            xbc = rtnorm(nsamples, mean[i], sd, a, b);
        } else if(lambda > 1e-6) {
            a = -1/lambda;
            b = infty;
            xbc = rtnorm(nsamples, mean[i], sd, a, b);            
        } else {
            xbc = rnorm(nsamples, mean[i], sd);
        }
        x = InvBoxCox(xbc, lambda);

        for(j=0; j<x.size(); j++) {
            EXz += fabs(x[j] - z[i]);
        }
        EXz = EXz/x.size();

        // Generate samples from Box Cox distribution
        if(lambda < -1e-6) {
            a = -infty;
            b = -1/lambda;
            xbc = rtnorm(nsamples, mean[i], sd, a, b);
        } else if(lambda > 1e-6) {
            a = -1/lambda;
            b = infty;
            xbc = rtnorm(nsamples, mean[i], sd, a, b);            
        } else {
            xbc = rnorm(nsamples, mean[i], sd);
        }
        x = InvBoxCox(xbc, lambda);

        nsamp = 0.0;
        for(j=0; j<(x.size()/2); j++) {
            EXXm += fabs(x[2*j] - x[2*j+1]);
            nsamp += 1.0;
        }
        EXXm = EXXm/nsamp;
        
        crps[i] = EXz - 0.5*EXXm;
    }
    return(crps);
}
