colnames(predictors.training) = paste("V", 1:dim(predictors.training)[2], sep = "")
colnames(predictors.test) = paste("V", 1:dim(predictors.test)[2], sep = "")

# Start with LASSO selection
cv = cv.glmnet(as.matrix(predictors.training), SWH.bc.standard.training, family = "gaussian", alpha = 1, nfold = 10)
lasso = glmnet(as.matrix(predictors.training), SWH.bc.standard.training, alpha = 1)
minindex = which.min(abs(lasso$lambda - cv$lambda.min))    
beta = lasso$beta[,minindex]
predictors.training = predictors.training[, which( abs(beta) > 1e-6 )]
predictors.test = predictors.test[, which( abs(beta) > 1e-6 )]

fit <- lm(SWH.bc.standard.training ~ ., data = predictors.training) # Use arima(...) later
fits = summary(fit)

# predict
SWH.bc.standard.pred <- predict(object = fit, newdata = predictors.test) #lm
SWH.bc.standard.pred.se <- fits$sigma

# Descale and detransform before computing different raknkings
SWH.bc.pred = SWH.bc.standard.pred*SWH.bc.sd.training + SWH.bc.mean.training
SWH.bc.pred.se = SWH.bc.standard.pred.se*SWH.bc.sd.training
SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambda.training) #detransform

SWH.bc.test = SWH.bc.standard.test*SWH.bc.sd.training + SWH.bc.mean.training
SWH.test = InvBoxCox(SWH.bc.test, SWH.bc.lambda.training) #detransfrom

# Handling NA's. Occurs rarely.
not.NA <- !is.na(SWH.test) & !is.na(SWH.pred)
if(sum(!not.NA ) > 0) {
    cat("j =", j, "k =", k, "sum(is.na(SWH.test)) =", sum(is.na(SWH.test)), "sum(is.na(SWH.pred)) =", sum(is.na(SWH.pred)), "\n")
}

U = ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambda.training)
save(U, file = paste("OutputComputations/Fourier_j_", j, "k_", k, ".Rdata", sep = ""))


idx.rank = idx.rank + 1
rankj[k,idx.rank,p,1] = mean(rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "MAE", list(distr = "boxcox", lambda = SWH.bc.lambda.training, se = SWH.bc.pred.se)), na.rm = TRUE)
rankj[k,idx.rank,p,2] = mean(rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "LOGS", list(distr = "boxcox", lambda = SWH.bc.lambda.training, se = SWH.bc.pred.se)), na.rm = TRUE)
rankj[k,idx.rank,p,3] = reliabilityIndex(ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambda.training), n.bins)
rankj[k,idx.rank,p,4] = reliabilityIndexSquare(ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambda.training), n.bins)
rankj[k,idx.rank,p,5] = RMSE(SWH.test[not.NA], SWH.bc.pred[not.NA], SWH.bc.pred.se, SWH.bc.lambda.training, n.bc.samp)
rankj[k,idx.rank,p,6] = mean(CRPS(SWH.test[not.NA], SWH.bc.pred[not.NA], SWH.bc.pred.se, SWH.bc.lambda.training, n.bc.samp), na.rm = TRUE)   
rankj[k,idx.rank,p,7] = mean(pSupport(mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambda.training), na.rm = TRUE)

#rankj.all[k,idx.rank,p,1,1:length(SWH.test[not.NA])] = rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "MAE", list(distr = "boxcox", lambda = SWH.bc.lambda.training, se = SWH.bc.pred.se))
#rankj.all[k,idx.rank,p,2,1:length(SWH.test[not.NA])] = rankModel(response = SWH.test[not.NA], fit, SWH.bc.pred[not.NA], model, "LOGS", list(distr = "boxcox", lambda = SWH.bc.lambda.training, se = SWH.bc.pred.se))
#rankj.all[k,idx.rank,p,3,1:length(SWH.test[not.NA])] = pnorm(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se))
#rankj.all[k,idx.rank,p,4,1:length(SWH.test[not.NA])] = CRPS(SWH.test[not.NA], SWH.bc.pred[not.NA], as.numeric(SWH.bc.pred.se), SWH.bc.lambda.training, n.bc.samp)

rankj.all[k,idx.rank,p,1,1:length(SWH.test[not.NA])] = ppredbc(SWH.bc.test[not.NA], mean = SWH.bc.pred[not.NA], sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambda.training)
rankj.all[k,idx.rank,p,2,1:length(SWH.test[not.NA])] = NA




#### Below when computing residuals. Note handle NA's a little different here then in the code above. Also note that we are only using the first part of the data set.
# predict
#SWH.bc.standard.pred <- predict(object = fit, newdata = predictors.training) #lm
#SWH.bc.standard.pred.se <- fits$sigma

# Descale and detransform before computing different raknkings
#SWH.bc.pred = SWH.bc.standard.pred*SWH.bc.sd.training + SWH.bc.mean.training
#SWH.bc.pred.se = SWH.bc.standard.pred.se*SWH.bc.sd.training
#SWH.pred = InvBoxCox(SWH.bc.pred, SWH.bc.lambda.training) #detransform

#SWH.bc.training = SWH.bc.standard.training*SWH.bc.sd.training + SWH.bc.mean.training
#SWH.training = InvBoxCox(SWH.bc.training, SWH.bc.lambda.training) #detransfrom

#bc.mean = meanboxcox(SWH.bc.pred, SWH.bc.pred.se, SWH.bc.lambda.training, n.bc.samp)
#rankj.all[k,idx.rank,p,1,idx.training] = SWH.training - bc.mean
#rankj.all[k,idx.rank,p,2,idx.training] = bc.mean
#rankj.all[k,idx.rank,p,3,idx.training] = qboxcox(0.5, mean = SWH.bc.pred, sd = SWH.bc.pred.se, lambda = SWH.bc.lambda.training)
#rankj.all[k,idx.rank,p,4,idx.training] = pnorm(SWH.bc.test, mean = SWH.bc.pred, sd = as.numeric(SWH.bc.pred.se))
#rankj.all[k,idx.rank,p,5,idx.training] = pSupport(mean = SWH.bc.pred, sd = SWH.bc.pred.se, lambda = SWH.bc.lambda.training)

#rankj.all[k,idx.rank,p,1,idx.training] = ppredbc(SWH.bc.test, mean = SWH.bc.pred, sd = as.numeric(SWH.bc.pred.se), lambda = SWH.bc.lambda.training)
#rankj.all[k,idx.rank,p,2,idx.training] = pSupport(mean = SWH.bc.pred, sd = SWH.bc.pred.se, lambda = SWH.bc.lambda.training)
