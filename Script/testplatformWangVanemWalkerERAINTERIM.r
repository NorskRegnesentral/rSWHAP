require(parallel)
require(glmnet)
require(forecast)
require(scoringRules)
require(Rcpp)
require(e1071)

sourceCpp("RMSE_CRPS_MEAN.cpp")


#options("width"=220)
#options("width"=180)

## Import functions to run the experiment below
source("functions.r")

n.training.test = 1

rank.methods = c("MAE", "LOGS", "RIDX", "RIDX2", "RMSE", "CRPS", "PSUP") 
#rank.methods = c("MAE", "LOGS", "RIDX", "RIDX2") 
n.rank.methods = length(rank.methods)

rank.methods.all = c("PIT", "PSUP") # Compute residuals
n.rank.methods.all = length(rank.methods.all)

area.longSWH = 1:181
area.latSWH = 1:68

#area.longSWH = 60:99
#area.latSWH = 22:45

n.bins = 20 # Number of bins when computing reliabillity index
n.bc.samp = 1e3 # Number of samples when estimating the Box Cox distribution (used to compute the reliability index)
bc.shift = 25 # To avoid negative values in the BC transform

seasons <- c("jfm", "amj", "jas", "ond")
n.seasons <- length(seasons)

mc.cores = min(length(area.longSWH), 32)

na.thresh = 500

run4seasons = FALSE
if(run4seasons) {
    # For the Wang models

    change.lambda = TRUE   # A new range is added to the BoxCoxLambda function
    
    models.4seasons = c("Wang et al. 2012", "Wang et al. 2012 no PC")
    n.models.4seasons = length(models.4seasons)

    n.time.max = 6000 # At least enough 
    rank.4seasons = array(dim = c(length(area.longSWH), length(area.latSWH), n.models.4seasons, n.training.test, n.rank.methods, n.seasons) )
    #rank.4seasons.all =array(dim = c(length(area.longSWH), length(area.latSWH), n.models.4seasons, n.training.test, 4, n.seasons, n.time.max) )

    # Wang et al. 2012 & Wang et al. 2012 LASSO & Wang no PC/EOF & Wang no PC/EOF seasonal
    for(s in 1:n.seasons) {
        load(paste("~/HDwaveData/Prepared/prepared_ERAInterimDaily_SWH_", seasons[s], ".RData", sep = "")) # On Research3
        ## Need to compute relative RMSE (relRMSE)
        SWH.baseline <- rowMeans(SWH[,,ind.baseline],dims = 2)

        # Compute different training and test parts
        training.test = list()
        part = list()

        part[[1]] = which(years %in% 1979:2000)
        part[[2]] = which(years %in% 2001:2015)
        training.test[[1]] = part

        # Denne returerer nå rank for hvert mål og n.latSWH. Denne vil jeg nå utvide til også å gå over alle modellene. Videre trenes alle modellene i tur i runParallellWangVanemWalkerModel
        rankj <- mclapply(area.longSWH, runParallelWangVanemWalker.4seasons, mc.cores = mc.cores, mc.preschedule=FALSE) 
        save(rankj, file = "~/HDwaveData/Results/rankjWangVanemWalkerERAINTERIM.4seasons.Rdata")
        #save(rankj.out, file = "rankj_out_4season.Rdata")
        for(j in 1:length(area.longSWH)) {
            rank.4seasons[j,,,,,s] = rankj[[j]][[1]][,,1,] # Use when n.training.test = 1
            #rank.4seasons.all[j,,,,,s,] = rankj[[j]][[2]][,,1,,] # Use when n.training.test = 1
        }
        
    }
    save(rank.4seasons, file = "OutputComputations/rankWangVanemWalkerERAINTERIM.4seaons.fullgrid.lambdamin_-0.5.Rdata")
    #save(rank.4seasons.all, file = "~/HDwaveData/Results/rankWangVanemWalkerERAINTERIM.all.4seaons.fullgrid.Rdata")
}

runFourier = TRUE
if(runFourier) {    
    # Data for the seasonal models
    load("~/HDwaveData/Prepared/SWH_SLP_SLP.grad.Rdata")
    load("Data/longlatSWHSLP.Rdata")

    spatial.neighborhoods = c(5) # For the Vanem&Walker model (max, min, mean).
    # models.fourier.bc below
    #models.fourier = c("SLP + SLP.grad", "SLP + SLP.grad + AR(2) (no Fourier)", "SLP + SLP.grad + AR(2)", "SLP + SLP.grad + neig = 5", "SLP + SLP.grad + AR(2) (no Fourier) neig = 5", "SLP + SLP.grad + AR(2) neig = 5")
    models.fourier = c("Proposed model, AR order = 2, spatial neigb. = 5", "Proposed model, AR order = 5, spatial neigb. = 5", "Proposed model, AR order = 10, spatial neigb. = 5")
    n.models.fourier = length(models.fourier)

    # For the Fourier transforms
    m.fourier = 365.25*4 # Number of observations per year. Four observations per day.
    K.fourier = 2 # Number of Fourier periodics
    n.time = length(time.all)

    intercept.fourier = fourier(rep(1,n.time), K=K.fourier, m=m.fourier)
    colnames(intercept.fourier) = paste("intercept", colnames(intercept.fourier), sep = "_")

    training.test = list()
    part = list()

    part[[1]] = which(years.all %in% 1979:2000)
    part[[2]] = which(years.all %in% 2001:2015)
    training.test[[1]] = part

    rank.fourier = array(dim = c(length(area.longSWH), length(area.latSWH), n.models.fourier, n.training.test, n.rank.methods) )
    #rank.fourier.all = array(dim = c(length(area.longSWH), length(area.latSWH), n.models.fourier, n.training.test, n.rank.methods.all, length(training.test[[1]][[2]])) )
    #rank.fourier.all = array(dim = c(length(area.longSWH), length(area.latSWH), n.models.fourier, n.training.test, n.rank.methods.all, length(training.test[[1]][[1]])) ) # When computing residuals

    rankj <- mclapply(area.longSWH, runParallelWangVanemWalker.fourier, mc.cores = mc.cores, mc.preschedule=FALSE)
    save(rankj, file = "~/HDwaveData/Results/rankjWangVanemWalkerERAINTERIM.fourier.fullgrid_Jarque–Bera.Rdata")
    #save(rankj, file = "~/HDwaveData/Results/residualsj.fourier.fullgrid.Rdata")
    for(j in 1:length(area.longSWH)) {
        cat("j =", j, "\n")
        rank.fourier[j,,,,] = rankj[[j]][[1]][,,1,] # Use when n.training.test = 1
        #rank.fourier.all[j,,,,,] = rankj[[j]][[2]][,,1,,] # Use when n.training.test = 1
    }
    save(rank.fourier, file = "OutputComputations/rankWangVanemWalkerERAINTERIM.fourier.fullgrid.lambdamin_0_Jarque–Bera.Rdata")    
    #save(rank.fourier.all, file = "~/HDwaveData/Results/rankWangVanemWalkerERAINTERIM.all.fourier.fullgrid.Rdata")
    #save(rank.fourier.all, file = "~/HDwaveData/Results/residuals.fourier.fullgrid.Rdata")
}
