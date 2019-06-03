##Validation of OOS curves

#This function returns a probability matrix for survival analysis models and 
#prediction vector for caret models
#If you want prediction vector of survival anlysis models run the function oos_predict
leave_one_out <- function(data, model, formula, weights = NULL, plot = TRUE, 
                          metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                          maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                          trControl = trainControl())
{
  len_d <- dim(data)[1]
  len_c <- length(unique(data$gap))
  if((model == "RandomSurvivalForest") || (model == "CoxPH"))
  {
    oosp <- matrix(0,nrow = len_d,ncol = len_c)
  }
  else{
    oosp <- c()
  }
  for(i in 1:len_d)
  {
    newd <<- data[-i,]
    nweight <<- weights[-i]
    if(model == "RandomSurvivalForest")
    {
      rfs <- rfsrc(formula,
                   data=newd, ntree = 100, samptype = "swr",
                   seed = 18,na.action = "na.impute", nodesize = 3, case.wt = nweight)
      newpred <- predict(rfs, newdata = data[i,],na.action = "na.impute")
      oosp[i,] <- c(1-newpred$survival,rep(1,(len_c-length(rfs$time.interest))))
      if(plot == TRUE)
      {
        if(i==1)
        {
          plot(rfs$time.interest,oosp[i,1:length(rfs$time.interest)], 'l', xlab = 'Time in days', 
               ylab = 'Out of stock probability', 
               main = 'Predicted out of stock probability curves RFS', ylim = c(0,1)) 
        }
        else
        {
          lines(rfs$time.interest,oosp[i,1:length(rfs$time.interest)])
        }
      }
    }
    else if(model == "CoxPH")
    {
      fit_cox <- coxph(formula, data = newd, weights = nweight)
      newpred <- survfit(fit_cox, newdata = data[i,])
      oosp[i,] <- c(1-newpred$surv,rep(1,(len_c-length(newpred$time))))
      if(plot == TRUE)
      {
        if(i==1)
        {
          plot(newpred$time,oosp[i,1:length(newpred$time)], 'l', xlab = 'Time in days', 
               ylab = 'Out of stock probability', 
               main = 'Predicted out of stock probability curves CoxPH', ylim = c(0,1)) 
        }
        else
        {
          lines(newpred$time,oosp[i,1:length(newpred$time)])
        }
      }
    }
    else
    {
      train_dim <- dim(newd)[2] - 2
      be <- train(x = newd[,1:train_dim], y = newd$gap, method = model, 
                  weights = nweight, metric = metric, maximize = maximize, 
                  trControl = trControl)
      newpred <- predict(be, newdata = data[i,])
      newpred <- ceiling(newpred)
      oosp <- c(oosp,newpred)
    }
  }
  if((model == "RandomSurvivalForest") || (model == "CoxPH"))
  {
    colnames(oosp) <- sort(unique(data$gap))
  }
  return(oosp)
}

#Converts probability matrix to predictions based on threshold
oospredict <- function(prob_matrix, threshold)
{
  pred <- c()
  len_d <- dim(prob_matrix)[1]
  for(i in 1:len_d)
  {
    dt <- prob_matrix[i,]
    if(anyNA(dt) == TRUE)
    {
      pred <- c(pred, NaN)
      #cat(NaN, " ")
      next
    }
    if(length(dt[dt >= threshold]) == 0)
    {
      ind1 <- as.integer(names(dt)[length(dt)])
      pred <- c(pred, ind1)
      next
    }
    ind1 <- as.integer(names(dt[dt >= threshold])[1]) #Day at which it is above threshold
    pr1 <- dt[dt >= threshold][1] #Probability at day ind1
    if(length(dt[dt < threshold]) == 0)
    {
      pred <- c(pred, ind1)
      #cat(ind1, " ")
      next
    }
    ind2 <- as.integer(names(dt[dt < threshold])[length(dt[dt < threshold])]) #Day at which it is below threshold
    pr2 <- dt[dt < threshold][length(dt[dt < threshold])] #Probability at day ind2
    
    #Using linear interpolation to find the day at which it crosses threshold
    day <- approx(x = c(pr1,pr2),y = c(ind1,ind2), xout = threshold)$y
    day <- ceiling(day)
    #cat(day, " ")
    pred <- c(pred,day)
  }
  return(pred)
}

rfday <- oospredict(prob_matrix = oosp, threshold = 0.8)
coxday <- oospredict(prob_matrix = oosp2, threshold = 0.8)


pr1 <- leave_one_out(data = model_ready_data, model = "RandomSurvivalForest", formula = Surv(gap, status) ~ 
                        promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                       +sales_rate_1+sales_rate_2+sales_rate_3)
exp1 <- oospredict(prob_matrix = pr1, threshold = 0.7)
pr2 <- leave_one_out(data = model_ready_data, model = "RandomSurvivalForest", formula = Surv(gap, status) ~ 
                       promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                     +sales_rate_1+sales_rate_2+sales_rate_3,weights = weight_train)
exp2 <- oospredict(prob_matrix = pr2, threshold = 0.7)
View(cbind(exp1,exp2,model_ready_data$gap))

cpr1 <- leave_one_out(data = model_ready_data, model = "CoxPH", formula = Surv(gap, status) ~ 
                       promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                      +sales_rate_1+sales_rate_2+sales_rate_3)
exp1 <- oospredict(prob_matrix = cpr1, threshold = 0.7)
cpr2 <- leave_one_out(data = model_ready_data, model = "CoxPH", formula = Surv(gap, status) ~ 
                       promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                      +sales_rate_1+sales_rate_2+sales_rate_3,weights = weight_train)
exp2 <- oospredict(prob_matrix = cpr2, threshold = 0.7)
View(cbind(exp1,exp2,model_ready_data$gap))

prob_bands <- function(prob_matrix, gap)
{
  len_d <- dim(prob_matrix)[1]
  prob_lr <- matrix(c(0,1), nrow = len_d, ncol = 2, byrow = TRUE)
  for(i in 1:len_d)
  {
    dt <- prob_matrix[i,]
    day <- gap[i]
    if(anyNA(dt) == TRUE)
    {
      prob_lr[i,] <- c(NaN, NaN)
      next
    }
    pr <- dt[as.character(day)]
    if(length(dt[dt < pr]) == 0)
    {
      pl<-0
      prob_lr[i,] <- c(pl,pr)
      next
    }
    prev_day <- as.integer(names(dt[dt < pr][length(dt[dt < pr])]))
    p_prev <- dt[dt < pr][length(dt[dt < pr])]
    if(prev_day == (day - 1))
    {
      pl <- p_prev
    }
    else
    {
      pl <- approx(x = c(prev_day,day),y = c(p_prev,pr), xout = (day-1))$y
    }
    prob_lr[i,] <- c(pl,pr)
  }
  return(prob_lr)
}

rfpband <- prob_bands(prob_matrix = oosp, gap = model_ready_data$gap)
coxpband <- prob_bands(prob_matrix = oosp2, gap = model_ready_data$gap)