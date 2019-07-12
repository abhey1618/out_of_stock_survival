##Validation of OOS curves

#This function returns a probability matrix for survival analysis models and 
#prediction vector for caret models
#If you want prediction vector of survival anlysis models run the function oos_predict
leave_one_out <- function(data, model, formula = NULL, formula2 = NULL, weights = NULL, plot = TRUE, 
                          metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                          maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                          trControl = trainControl(), tau = 0.5)
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
               main = 'Predicted out of stock probability curves RFS', ylim = c(0,1),xlim = c(0,50)) 
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
               main = 'Predicted out of stock probability curves CoxPH', ylim = c(0,1), xlim = c(0,50)) 
        }
        else
        {
          lines(newpred$time,oosp[i,1:length(newpred$time)])
        }
      }
    }
    else if(model == "QuantReg")
    {
      rqfit <- rq(formula = formula2, data = newd, tau = tau)
      newpred <- predict(rqfit, newdata = data[i,])
      oosp <- c(oosp, newpred)
    }
    else
    {
      train_dim <- dim(newd)[2] - 3
      be <- train(form = formula2, data = newd, method = model, 
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

leave_one_out_cnsr <- function(data, model, formula, weights = NULL, plot = TRUE, 
                          metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                          maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                          trControl = trainControl(), tau = 0.5)
{
  len_d <- dim(data)[1]
  len_r <- sum(data$status == 1)
  len_c <- length(unique(data$gap))
  if((model == "RandomSurvivalForest") || (model == "CoxPH"))
  {
    oosp <- matrix(0,nrow = len_r,ncol = len_c)
  }
  else{
    oosp <- c()
  }
  plt <- 0
  j<-0
  for(i in 1:len_d)
  {
    # For first plot 
    if(data$status[i] == 0) next
    j <- j+1
    newd <<- data[-i,]
    nweight <<- weights[-i]
    if(model == "RandomSurvivalForest")
    {
      rfs <- rfsrc(formula,data=newd, ntree = 100, samptype = "swr",
                   seed = 18,na.action = "na.impute", nodesize = 3, case.wt = nweight)
      newpred <- predict(rfs, newdata = data[i,],na.action = "na.impute")
      oosp[j,] <- c(1-newpred$survival,rep(1,(len_c-length(rfs$time.interest))))
      if(plot == TRUE)
      {
        if(plt==0)
        {
          plt = 1
          plot(rfs$time.interest,oosp[j,1:length(rfs$time.interest)], 'l', xlab = 'Time in days', 
               ylab = 'Out of stock probability', 
               main = 'Predicted out of stock probability curves RFS', ylim = c(0,1)) 
        }
        else
        {
          lines(rfs$time.interest,oosp[j,(1:length(rfs$time.interest))])
        }
      }
    }
    else if(model == "CoxPH")
    {
      fit_cox <- coxph(formula, data = newd, weights = nweight)
      newpred <- survfit(fit_cox, newdata = data[i,])
      oosp[j,] <- c(1-newpred$surv,rep(1,(len_c-length(newpred$time))))
      if(plot == TRUE)
      {
        if(plt==0)
        {
          plt <- 1
          plot(newpred$time,oosp[j,1:length(newpred$time)], 'l', xlab = 'Time in days', 
               ylab = 'Out of stock probability', 
               main = 'Predicted out of stock probability curves CoxPH', ylim = c(0,1), 
               xlim = c(0,max(newpred$time))) 
        }
        else
        {
          lines(newpred$time,oosp[j,1:length(newpred$time)])
        }
      }
    }
    else if(model == "QuantReg")
    {
      rqfit <- rq(formula = formula, data = newd, tau = tau)
      newpred <- predict(rqfit, newdata = data[i,])
      oosp <- c(oosp, newpred)
    }
    else
    {
      train_dim <- dim(newd)[2] - 3
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

#Function to create a data frame for insample predictions for the 10 models

#Function to create insample predictions weighted models
create_insample_predictions_df_w <- function(formula, formula2, models, data, weight,
                                             metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                                             maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                                             trControl = trainControl())
{
  wei <<- weight
  data_inp <<- data
  df <- data.frame(observed = data$gap)
  len_m <- length(models[[1]])
  ##cat("len"," " ,len_m, "\n")
  for(i in 1:len_m)
  {
    #cat(i, "\n")
    #ifelse(anyNA(models[[2]][,i]), wei <- NULL,wei <- models[[2]][,i]) 
    model = models[[1]][i]
    if(model == "RandomSurvivalForest")
    {
      #cat("Yeah rfs")
      rfs <- rfsrc(formula = formula, data=data_inp, ntree = 100, 
                   samptype = "swr",seed = 18, nodesize = 3, case.wt = weight)
      pr1 <- 1 - rfs$survival
      colnames(pr1) <- rfs$time.interest
      p <- models[[2]][i]
      df[paste("RFS",p,"w", sep = "_")] <- oospredict(pr1, p)
    }
    else if(model == "CoxPH")
    {
      #cat("Yeah CoxPH")
      cox_fit <- coxph(formula = formula, data=data_inp, weights = wei)
      pr1 <- t(1 - survfit(cox_fit, newdata = data_inp)$surv)
      colnames(pr1) <- survfit(cox_fit, newdata = data_inp)$time
      p <- models[[2]][i]
      df[paste("CoxPH",p,"w", sep = "_")] <- oospredict(pr1, p)
    }
    else
    {
      train_dim <- dim(data)[2] - 3
      #cat("Yeah Caret")
      newm<- train(form = formula2, data = data_inp, method = model, 
                   weights = wei, metric = metric, maximize = maximize, 
                   trControl = trControl)
      newpred <- predict(newm, newdata = data[,1:train_dim])
      newpred <- ceiling(newpred)
      df[paste(model,"w",sep = "_")] <- newpred
    }
  }
  return(df)
}

#Creates in sample predictions for unweighted models
create_insample_predictions_df <- function(formula, formula2, models, data,
                                           metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                                           maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                                           trControl = trainControl())
{
  data_inp <<- data
  df <- data.frame(observed = data$gap)
  len_m <- length(models[[1]])
  #cat("len"," " ,len_m, "\n")
  for(i in 1:len_m)
  {
    #cat(i, "\n")
    #ifelse(anyNA(models[[2]][,i]), wei <- NULL,wei <- models[[2]][,i]) 
    model = models[[1]][i]
    if(model == "RandomSurvivalForest")
    {
      #cat("Yeah rfs")
      rfs <- rfsrc(formula = formula, data=data, ntree = 100, 
                   samptype = "swr",seed = 18, nodesize = 3)
      pr1 <- 1 - rfs$survival
      colnames(pr1) <- rfs$time.interest
      p <- models[[2]][i]
      df[paste("RFS",p, sep = "_")] <- oospredict(pr1, p)
    }
    else if(model == "CoxPH")
    {
      #cat("Yeah CoxPH")
      cox_fit <- coxph(formula = formula,data=data_inp)
      pr1 <- t(1 - survfit(cox_fit, newdata = data_inp)$surv)
      colnames(pr1) <- survfit(cox_fit, newdata = data_inp)$time
      p <- models[[2]][i]
      df[paste("CoxPH",p, sep = "_")] <- oospredict(pr1, p)
    }
    else
    {
      train_dim <- dim(data)[2] - 3
      #cat("Yeah Caret")
      newm<- train(form = formula2, data = data_inp, method = model, 
                   metric = metric, maximize = maximize, trControl = trControl)
      newpred <- predict(newm, newdata = data[,1:train_dim])
      newpred <- ceiling(newpred)
      df[paste(model)] <- newpred
    }
  }
  return(df)
}

#Function to create out of sample predictions old models
create_oosample_predictions_df_w_old <- function(formula, formula2,models, data, weights,
                                             metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                                             maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                                             trControl = trainControl())
{
  df <- data.frame(observed = data$gap)
  len_m <- length(models[[1]])
  cat("len"," " ,len_m, "\n")
  for(i in 1:len_m)
  {
    cat(i, "\n")
    #ifelse(anyNA(models[[2]][,i]), wei <- NULL,wei <- models[[2]][,i]) 
    model = models[[1]][i]
    if(model == "RandomSurvivalForest")
    {
      cat("Yeah rfs")
      pr1 <- leave_one_out(data = data, model = model, formula = formula, weights = weights, 
                           plot = FALSE)
      p <- models[[2]][i]
      df[paste("RFS",p,"w")] <- oospredict(pr1, p)
    }
    else if(model == "CoxPH")
    {
      cat("Yeah CoxPH")
      pr1 <- leave_one_out(data = data, model = model, formula = formula, weights = weights, 
                           plot = FALSE)
      p <- models[[2]][i]
      df[paste("CoxPH",p,"w")] <- oospredict(pr1, p)
    }
    else
    {
      cat("Yeah Caret")
      newpred <- leave_one_out(data = data, model = model, formula2 = formula2, weights = weights, 
                               plot = FALSE, metric = metric, maximize = maximize, trControl = trControl)
      #newpred <- ceiling(newpred)
      df[paste(model,"w")] <- newpred
    }
  }
  return(df)
}

#Creates out of sample leave one out predictions for unweighted models
create_oosample_predictions_df_old <- function(formula,formula2, models, data,
                                           metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                                           maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                                           trControl = trainControl())
{
  df <- data.frame(observed = data$gap)
  len_m <- length(models[[1]])
  cat("len"," " ,len_m, "\n")
  for(i in 1:len_m)
  {
    cat(i, "\n")
    #ifelse(anyNA(models[[2]][,i]), wei <- NULL,wei <- models[[2]][,i]) 
    model = models[[1]][i]
    if(model == "RandomSurvivalForest")
    {
      cat("Yeah rfs")
      pr1 <- leave_one_out(data = data, model = model, formula = formula, 
                           plot = FALSE)
      p <- models[[2]][i]
      df[paste("RFS",p)] <- oospredict(pr1, p)
    }
    else if(model == "CoxPH")
    {
      cat("Yeah CoxPH")
      pr1 <- leave_one_out(data = data, model = model, formula = formula, 
                           plot = FALSE)
      p <- models[[2]][i]
      df[paste("CoxPH",p)] <- oospredict(pr1, p)
    }
    else
    {
      cat("Yeah Caret")
      newpred <- leave_one_out(data = data, model = model, formula2 = formula2,
                               plot = FALSE, metric = metric, maximize = maximize, trControl = trControl)
      #newpred <- ceiling(newpred)
      df[paste(model)] <- newpred
    }
  }
  return(df)
}

create_oosample_predictions_df_w <- function(formula, formula2,models, data, weights,
                                             metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                                             maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                                             trControl = trainControl())
{
  df <- data.frame(observed = data$gap[data$status == 1])
  len_m <- length(models[[1]])
  #cat("len"," " ,len_m, "\n")
  for(i in 1:len_m)
  {
    #cat(i, "\n")
    #ifelse(anyNA(models[[2]][,i]), wei <- NULL,wei <- models[[2]][,i]) 
    model = models[[1]][i]
    if(model == "RandomSurvivalForest")
    {
      #cat("Yeah rfs")
      pr1 <- leave_one_out_cnsr(data = data, model = model, formula = formula, weights = weights, 
                                plot = FALSE)
      p <- models[[2]][i]
      df[paste("RFS",p,"w", sep = "_")] <- oospredict(pr1, p)
    }
    else if(model == "CoxPH")
    {
      #cat("Yeah CoxPH")
      pr1 <- leave_one_out_cnsr(data = data, model = model, formula = formula, weights = weights, 
                                plot = FALSE)
      p <- models[[2]][i]
      df[paste("CoxPH",p,"w", sep = "_")] <- oospredict(pr1, p)
    }
    else
    {
      #cat("Yeah Caret")
      newpred <- leave_one_out(data = data, model = model, formula2 = formula2, weights = weights, 
                               plot = FALSE, metric = metric, maximize = maximize, trControl = trControl)
      #newpred <- ceiling(newpred)
      df[paste(model,"w")] <- newpred
    }
  }
  return(df)
}

create_oosample_predictions_df <- function(formula,formula2, models, data,
                                           metric = ifelse(is.factor(y), "Accuracy", "RMSE"),
                                           maximize = ifelse(metric == "RMSE", FALSE, TRUE),
                                           trControl = trainControl())
{
  df <- data.frame(observed = data$gap[data$status == 1])
  len_m <- length(models[[1]])
  #cat("len"," " ,len_m, "\n")
  for(i in 1:len_m)
  {
    #cat(i, "\n")
    #ifelse(anyNA(models[[2]][,i]), wei <- NULL,wei <- models[[2]][,i]) 
    model = models[[1]][i]
    if(model == "RandomSurvivalForest")
    {
      #cat("Yeah rfs")
      pr1 <- leave_one_out_cnsr(data = data, model = model, formula = formula, 
                                plot = FALSE)
      p <- models[[2]][i]
      df[paste("RFS",p, sep = "_")] <- oospredict(pr1, p)
    }
    else if(model == "CoxPH")
    {
      #cat("Yeah CoxPH")
      pr1 <- leave_one_out_cnsr(data = data, model = model, formula = formula, 
                                plot = FALSE)
      p <- models[[2]][i]
      df[paste("CoxPH",p, sep = "_")] <- oospredict(pr1, p)
    }
    else
    {
      #cat("Yeah Caret")
      newpred <- leave_one_out(data = data, model = model, formula2 = formula2,
                               plot = FALSE, metric = metric, maximize = maximize, trControl = trControl)
      #newpred <- ceiling(newpred)
      df[paste(model)] <- newpred
    }
  }
  return(df)
}
