##Validation of OOS curves

leave_one_out <- function(data, model)
{
  len_d <- dim(data)[1]
  len_c <- length(unique(data$gap))
  oosp <- matrix(0,nrow = len_d,ncol = len_c)
  for(i in 1:len_d)
  {
    newd <- data[-i,]
    if(model == "RandomSurvivalForest")
    {
      rfs <- rfsrc(Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2,
                   data=newd, ntree = 100, samptype = "swr",
                   seed = 18,na.action = "na.impute", nodesize = 5)
      newpred <- predict(rfs, newdata = data[i,],na.action = "na.impute")
      oosp[i,] <- c(1-newpred$survival,rep(1,(len_c-length(rfs$time.interest))))
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
    else if(model == "CoxPH")
    {
      fit_cox <- coxph(Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2,
                   data=newd)
      newpred <- survfit(fit_cox, newdata = data[i,])
      oosp[i,] <- c(1-newpred$surv,rep(1,(len_c-length(newpred$time))))
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
  colnames(oosp) <- sort(unique(data$gap))
  return(oosp)
}

oosp <- leave_one_out(data = model_ready_data, model = "RandomSurvivalForest")
oosp2 <- leave_one_out(data = model_ready_data, model = "CoxPH")

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