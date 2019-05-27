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
             ylab = 'Out of stock probability', main = 'Predicted out of stock probability curves', ylim = c(0,1)) 
      }
      else
      {
        lines(rfs$time.interest,oosp[i,1:length(rfs$time.interest)])
      }
    }
  }
  colnames(oosp) <- sort(unique(data$gap))
  return(oosp)
}

oosp <- leave_one_out(data = model_ready_data, model = "RandomSurvivalForest")
