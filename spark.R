library(data.table)
library(MTS)
library(sqldf)
library(dplyr)
library(survival)
library(ggplot2)
#library(ggfortify)
#library(randomForestSRC)
#library(coxme)
library(quantreg)
library(SparkR)
library(caret)
library(chemometrics)

schema_rdd <- "sku INT, p_sls_d STRING, p_week_start_date STRING, co_loc_i INT, 
istk_appl_f INT, sum_elig_istk_out_of_stk_f INT, circ_f DOUBLE, tpc_f DOUBLE, 
baseprice DOUBLE, offer_price_orig DOUBLE, exit_price DOUBLE, offer_price DOUBLE,
nrgtn INT, ttl_units INT, ttl_units_gross INT, dy_itm_loc_oos_ind INT, sum_boh_q DOUBLE,
sum_eoh_q DOUBLE, promo_ind INT, item_loc STRING"

all_itm_oos_agg <- read.df("./oos_sample_data_2.csv",source = "csv",header=T,
                           schema = schema_rdd)
# all_itm_oos_agg_df <- fread("./oos_sample_data_2.csv",header=T, stringsAsFactors = F)
# unique_itm_loc <- unique(all_itm_oos_agg_df$item_loc)
# itm_loc_smpl <- unique_itm_loc[1:100]
# all_itm_oos_agg_smpl <- all_itm_oos_agg_df[all_itm_oos_agg_df$item_loc%in%itm_loc_smpl,]
# 
# all_itm_oos_agg <- as.DataFrame(all_itm_oos_agg_smpl)

str(all_itm_oos_agg)

newsparkfunc <- function(key,x)
{
  library(data.table)
  library(MTS, lib.loc = "/home_dir/z004189/Library/R")
  library(sqldf)
  library(dplyr)
  library(survival)
  library(quantreg)
  library(SparkR)
  library(caret)
  library(chemometrics, lib.loc = "/home_dir/z004189/Library/R")
  
  #Locate runs
  run_locate <- function(y){
    start_vec <- c()
    end_vec <- c()
    for(i in 2:length(y)){
      if(y[i]==1 & y[i-1]==0) start_vec <- c(start_vec,i)
      if((y[i-1]==1 & i==2)) start_vec <- c(start_vec,i-1)
      
      if(y[i]==0 & y[i-1]==1 & i <= length(y)) end_vec <- c(end_vec,(i-1))
      if(y[i]==1 & i == length(y)) end_vec <- c(end_vec,i)
    }
    return(data.frame(start_vec=start_vec,end_vec=end_vec))
  }
  
  #Locate increase in inventories
  inv_up_locate <- function(z,y){
    start_vec <- c()
    end_vec <- c()
    mark1 <- 0
    mark2 <- 0
    for(i in 2:length(z)){
      #mark1 is to remember start of an oos
      #mark2 is to see if an oos is still going on or not
      if(y[i] == 1 & mark2 == 0)
      {
        mark1 <- i
        mark2 <- 1
      }
      if((z[i]>0 & i==2)) start_vec <- c(start_vec,1)
      if(z[i] > z[i-1] & y[i] == 0)
      {
        start_vec <- c(start_vec,i)
        if(y[i-1] == 0)
        {
          end_vec <- c(end_vec,(i-1))
        }
        else
        {
          end_vec <- c(end_vec, mark1)
          mark2 <- 0
        }
      }
      
      #if(z[i]==0 & z[i-1]==1 & i <= length(z)) end_vec <- c(end_vec,(i-1))
      #if(z[i]==1 & i == length(z)) end_vec <- c(end_vec,i)
    }
    end_vec <- c(end_vec,length(z))
    if(end_vec[1] == 0)
    {
      start_vec <- start_vec[-1]
      end_vec <- end_vec[-1]
    }
    if(end_vec[1] < start_vec[1])
    {
      end_vec <- end_vec[-1]
    }
    return(data.frame(start_vec=start_vec,end_vec=end_vec))
  }
  
  #Function to generate covariates
  promo_extractor <- function(inp_dt,lg){
    #inp_dt <- all_itm_oos_agg_949663_loc_aggrgtd[1:5,]
    #lg <- 28
    ln_x = dim(inp_dt)[1]
    int <- as.integer(lg/4) #Interval length of each of the 4 parts
    oos_smry_1 <- rep(0,ln_x)
    oos_smry_2 <- rep(0,ln_x)
    promo_smry_1 <- rep(0,ln_x)
    promo_smry_2 <- rep(0,ln_x)
    sales_smry <- rep(0, ln_x)
    sales_rate_1 <- rep(0,ln_x)
    sales_rate_2 <- rep(0, ln_x)
    sales_rate_3 <- rep(0, ln_x)
    #refill_1 <- rep(0,ln_x)
    #refill_2 <- rep(0,ln_x)
    for(i in 1:ln_x){
      ##cat(" ###   I   ###" , i, "\n")
      
      if(i <= lg){
        z <- inp_dt[(1:ifelse(i==1,1,i-1)),]
        #z$promo_ind <- as.integer(z$promo_ind)
        promo_smry_1[i] <- ifelse(sum(z$promo_ind)==0,0,(sum(z$promo_ind)/min(lg,i)))
        promo_smry_2[i] <- mean(z$baseprice[z$promo_ind==1] - z$offer_price[z$promo_ind==1])
      }
      if(i > lg){
        z <- inp_dt[(i-lg):(i-1),]   # assuming lg >= 2
        #z$promo_ind <- as.integer(z$promo_ind)
        promo_smry_1[i] <- ifelse(sum(z$promo_ind)==0,0,(sum(z$promo_ind)/min(lg,i)))
        promo_smry_2[i] <- mean(z$baseprice[z$promo_ind==1] - z$offer_price[z$promo_ind==1])
      }
      
      # if(i > aa$start_vec[1]){
      
      #
      #
      oos_smry_1[i] <- ifelse(sum(z$dy_itm_loc_oos_ind)==0,0,(sum(z$dy_itm_loc_oos_ind)/min(lg,i)))
      sales_smry[i] <- sum(z$ttl_units)
      
      #Dividing the df z into 4 parts to see the rate of change of sales over weeks
      if(i <= int)
      {
        sales_rate_1[i] <- 0
        sales_rate_2[i] <- 0
        sales_rate_3[i] <- 0
      }
      else if(i <= (2*int))
      {
        s1 <- sum(z$ttl_units[1:int])
        s2 <- sum(z$ttl_units) - s1
        sales_rate_1[i] <- (s2-s1)/s1
        sales_rate_2[i] <- 0
        sales_rate_3[i] <- 0
      }
      else if(i <= (3*int))
      {
        s1 <- sum(z$ttl_units[1:int])
        s2 <- sum(z$ttl_units[1:(2*int)]) - s1
        s3 <- sum(z$ttl_units) - s1 - s2
        sales_rate_1[i] <- (s2-s1)/s1
        sales_rate_2[i] <- (s3-s2)/s2
        sales_rate_3[i] <- 0
      }
      else
      {
        s1 <- sum(z$ttl_units[1:int])
        s2 <- sum(z$ttl_units[1:(2*int)]) - s1
        s3 <- sum(z$ttl_units[1:(3*int)]) - s1 - s2
        s4 <- sum(z$ttl_units) - s1 - s2 - s3
        sales_rate_1[i] <- (s2-s1)/s1
        sales_rate_2[i] <- (s3-s2)/s2
        sales_rate_3[i] <- (s4-s3)/s3
      }
      
      y <- z$dy_itm_loc_oos_ind
      start_vec <- c()
      end_vec <- c()
      if(length(y)==1) {start_vec <- c(start_vec,1) ; end_vec <- c(end_vec,1)}
      if(length(y) > 1){
        for(j in 2:length(y)){
          if(y[j]==1 & y[j-1]==0) start_vec <- c(start_vec,j)
          if((y[j]==1 & y[j-1]==1 & j==2)) start_vec <- c(start_vec,(j-1))
          if((y[j]==0 & y[j-1]==1 & j==2)) start_vec <- c(start_vec,(j-1))
          
          if(y[j]==0 & y[j-1]==1 & j <= length(y)) end_vec <- c(end_vec,(j-1))
          if(y[j]==1 & y[j-1]==0 & j == length(y)) end_vec <- c(end_vec,j)
          if(y[j]==1 & y[j-1]==1 & j == length(y)) end_vec <- c(end_vec,j)
          
        }
      }
      oos_smry_2[i] <-  ifelse(length(end_vec) ==0 & length(start_vec)==0,0,max(end_vec - start_vec)+1)  
    }
    return(data.frame(oos_smry_1=oos_smry_1,oos_smry_2=oos_smry_2,promo_smry_1=promo_smry_1,
                      promo_smry_2=promo_smry_2, sales_smry = sales_smry, sales_rate_1 = sales_rate_1,
                      sales_rate_2 = sales_rate_2, sales_rate_3 = sales_rate_3))
  }
  
  #Create data for model with censoring
  create_data_for_model_new <- function(dataframe, lg)
  {
    aa <- inv_up_locate(dataframe$sum_boh_q,dataframe$dy_itm_loc_oos_ind)
    #aa$start_0 <- lag(aa$end_vec,1)
    #aa$start_0[1] <- 0
    bb <- promo_extractor(dataframe,lg)
    bb<-cbind(bb, dataframe$promo_ind)
    colbb <- dim(bb)[2]
    colnames(bb)[colbb] <- "promo_ind"
    bb<-cbind(bb, dataframe$sum_boh_q)
    colbb <- dim(bb)[2]
    colnames(bb)[colbb] <- "sum_boh_q"
    
    #abb <- cbind(bb,lag(dataframe$dy_itm_loc_oos_ind,1))
    model_ready_data <- bb[aa$start_vec,]
    model_ready_data$status <- dataframe$dy_itm_loc_oos_ind[aa$end_vec]
    model_ready_data$gap <- aa$end_vec - aa$start_vec + 1
    # new<-cbind(tail(abb,1), (dim(abb)[1] - tail(aa$end_vec,1)))
    # colnames(new)[(colbb+2)] <- "gap"
    # model_ready_data <- rbind(model_ready_data, new)
    model_ready_data <- model_ready_data[-1,]
    model_ready_data$sku <- rep(dataframe$sku[1], dim(model_ready_data)[1])
    
    return(model_ready_data)
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
  
  dist_to_weight <- function(distance)
  {
    weight <- 1/(distance+1)
    weight <- weight/sum(weight)
    return(weight)
  }
  
  appl_weight <- function(differ)
  {
    ifelse((differ > 0), 1.2, 0.8)
  }
  
  my_metric <- function(data, lev = NULL, model = NULL) 
  {
    diff <- data$pred - data$obs
    met <- sum(((abs(diff)/data$obs) * sapply(diff,appl_weight)))
    c(MYM = met)
  }
  
  wtd_metric <- function(expected, observed)
  {
    diff <- expected - observed
    cat("Exp",expected, "\n")
    cat("Obs",observed, "\n")
    cat("NA in expected ->", anyNA(expected), "\n")
    cat("NA in observed ->", anyNA(observed), "\n")
    cat("0 in expected ->", sum(expected == 0), "\n")
    met <- mean(((abs(diff)/expected) * sapply(diff,appl_weight)))
    return(met)
  }
  
  calc_metric <- function(df)
  {
    df <- as.data.frame(df)
    metric <- c()
    obs <- df[,1]
    for(i in 2:dim(df)[2])
    {
      met <- wtd_metric(expected = df[,i], observed = obs)
      metric <- c(metric,met)
      names(metric)[(i-1)] <- colnames(df)[i]
    }
    cat("calc_metric", metric, "\n")
    return(metric)
  }
  
  cat("START", "\n")
  
  # start_vec <- c()
  # end_vec <- c()
  # inp_x <- x[,"dy_itm_loc_oos_ind"]
  # inp_x <- as.data.frame(inp_x)
  # for(i in 2:length(inp_x)){
  #   if(inp_x[i]==1 & inp_x[i-1]==0) start_vec <- c(start_vec,i)
  #   if((inp_x[i-1]==1 & i==2)) start_vec <- c(start_vec,i-1)
  #   
  #   if(inp_x[i]==0 & inp_x[i-1]==1 & i <= length(inp_x)) end_vec <- c(end_vec,(i-1))
  #   if(inp_x[i]==1 & i == length(inp_x)) end_vec <- c(end_vec,i)
  # }
  # 
  # model_ready_data <- data.frame(start_vec=start_vec,end_vec=end_vec)
  
  #print(dim(model_ready_data))
  #cat(class(x), "\n")
  #x <- SparkR:::collect(x)
  #cat("END", "\n")
  
  model_ready_data <- try(create_data_for_model_new(x, 30));
  if(class(model_ready_data) == "try-error") model_ready_data <- NULL
  #if(dim(model_ready_data)[1] == 0) model_ready_data <- NULL
  if(is.null(model_ready_data))
  {
    # metos <- rep(100, 4)
    # names(metos) <- c("CoxPH_0.4_w","CoxPH_0.3_w","CoxPH_0.4","CoxPH_0.3")
    metos <- rep(100, 2)
    names(metos) <- c("CoxPH_0.4","CoxPH_0.3")
  }
  else if(sum(model_ready_data$status) > 5)
  {
    model_ready_data <- na.omit(model_ready_data)
    cat("NA's ommited", "\n")
    
    train_dim <- dim(model_ready_data)[2] - 3
    
    #Robust Mahalanobis Distance for weighting
    #Doing this over 1:(train_dim - 1) because can ignore promo_ind for weighting
    # weight_train <- try(dist_to_weight(Moutlier(model_ready_data[,2:(train_dim-2)]
    #                                                 ,plot = FALSE)$rd))
    # 
    # if(class(weight_train) == "try-error") weight_train <- rep(1,dim(model_ready_data)[1])
    # cat("Assigned weights", "\n")
    
    thresh <- c(0.4,0.3)
    l <- list("lmodels" = c("CoxPH","CoxPH"), 
              "thresholds" = rep(thresh,1))
    # looval_pred <- create_oosample_predictions_df_w(formula = Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+
    #                                               promo_smry_2+sales_smry+sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q, models = l,
    #                                             data = model_ready_data, weight = weight_train,
    #                                             metric = "MYM", maximize = FALSE, 
    #                                             trControl = trainControl(summaryFunction = my_metric))
    cat(x$item_loc[1], "\n")
    cat('Weighted models trained', "\n")
    # l <- list("lmodels" = c("CoxPH","CoxPH","RandomSurvivalForest","RandomSurvivalForest"), 
    #           "thresholds" = rep(thresh,2))
    
    looval_pred2 <- try(create_oosample_predictions_df(formula = Surv(gap, status) ~ promo_ind+oos_smry_1+
                                                     oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                                                   +sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q,
                                                   models = l,data = model_ready_data,
                                                   metric = "MYM", maximize = FALSE, 
                                                   trControl = trainControl(summaryFunction = my_metric)))
    
    cat('Unweighted models trained', "\n")
    #looval_pred <- cbind(looval_pred,looval_pred2[-1])
    cat('Results combined', "\n")
    if(class(looval_pred2) == "try-error")
    {
      metos <- rep(150, 2)
      names(metos) <- c("CoxPH_0.4","CoxPH_0.3")
    }
    else metos <- calc_metric(df = looval_pred2)
    cat('metos ready', "\n")
    #res <- data.frame(key, dim(model_ready_data)[1], stringsAsFactors = FALSE)
    #colnames(res) <- c("key", "numrows")
  }
  else if(sum(model_ready_data$status) > 1)
  {
    model_ready_data <- na.omit(model_ready_data)
    cat("Else 2 case", x$item_loc[1], "\n")
    looval <- data.frame(observed = model_ready_data$gap[model_ready_data$status == 1])
    # looval$CoxPH_0.4_w <- rep(mean(looval$observed),dim(looval)[1])
    # looval$CoxPH_0.3_w <- rep(mean(looval$observed),dim(looval)[1])
    looval$CoxPH_0.4 <- rep(mean(looval$observed),dim(looval)[1])
    looval$CoxPH_0.3 <- rep(mean(looval$observed),dim(looval)[1])
    
    metos <- calc_metric(df = looval)
  }
  else
  {
    metos <- rep(100, 2)
    cat("0 oos case", x$item_loc[1], "\n")
    names(metos) <- c("CoxPH_0.4","CoxPH_0.3")
  }
  #z <- data.frame(itm_loc=key,data_row=dim(y)[1], stringsAsFactors = FALSE)
  #return(y)
  res <- data.frame(t(c(key, metos)), stringsAsFactors = FALSE)
  colnames(res)[1] <- "Key"
  res
}

result <- gapplyCollect(all_itm_oos_agg, cols = all_itm_oos_agg$item_loc, newsparkfunc)

newsparkfunc2 <- function(key,x)
{
  library(data.table)
  library(MTS, lib.loc = "/home_dir/z004189/Library/R")
  library(sqldf)
  library(dplyr)
  library(survival)
  library(quantreg)
  library(SparkR)
  library(caret)
  library(chemometrics, lib.loc = "/home_dir/z004189/Library/R")
  
  #Locate runs
  run_locate <- function(y){
    start_vec <- c()
    end_vec <- c()
    for(i in 2:length(y)){
      if(y[i]==1 & y[i-1]==0) start_vec <- c(start_vec,i)
      if((y[i-1]==1 & i==2)) start_vec <- c(start_vec,i-1)
      
      if(y[i]==0 & y[i-1]==1 & i <= length(y)) end_vec <- c(end_vec,(i-1))
      if(y[i]==1 & i == length(y)) end_vec <- c(end_vec,i)
    }
    return(data.frame(start_vec=start_vec,end_vec=end_vec))
  }
  
  #Locate increase in inventories
  inv_up_locate <- function(z,y){
    start_vec <- c()
    end_vec <- c()
    mark1 <- 0
    mark2 <- 0
    for(i in 2:length(z)){
      #mark1 is to remember start of an oos
      #mark2 is to see if an oos is still going on or not
      if(y[i] == 1 & mark2 == 0)
      {
        mark1 <- i
        mark2 <- 1
      }
      if((z[i]>0 & i==2)) start_vec <- c(start_vec,1)
      if(z[i] > z[i-1] & y[i] == 0)
      {
        start_vec <- c(start_vec,i)
        if(y[i-1] == 0)
        {
          end_vec <- c(end_vec,(i-1))
        }
        else
        {
          end_vec <- c(end_vec, mark1)
          mark2 <- 0
        }
      }
      
      #if(z[i]==0 & z[i-1]==1 & i <= length(z)) end_vec <- c(end_vec,(i-1))
      #if(z[i]==1 & i == length(z)) end_vec <- c(end_vec,i)
    }
    end_vec <- c(end_vec,length(z))
    if(end_vec[1] == 0)
    {
      start_vec <- start_vec[-1]
      end_vec <- end_vec[-1]
    }
    if(end_vec[1] < start_vec[1])
    {
      end_vec <- end_vec[-1]
    }
    return(data.frame(start_vec=start_vec,end_vec=end_vec))
  }
  
  #Function to generate covariates
  promo_extractor <- function(inp_dt,lg){
    #inp_dt <- all_itm_oos_agg_949663_loc_aggrgtd[1:5,]
    #lg <- 28
    ln_x = dim(inp_dt)[1]
    int <- as.integer(lg/4) #Interval length of each of the 4 parts
    oos_smry_1 <- rep(0,ln_x)
    oos_smry_2 <- rep(0,ln_x)
    promo_smry_1 <- rep(0,ln_x)
    promo_smry_2 <- rep(0,ln_x)
    sales_smry <- rep(0, ln_x)
    sales_rate_1 <- rep(0,ln_x)
    sales_rate_2 <- rep(0, ln_x)
    sales_rate_3 <- rep(0, ln_x)
    #refill_1 <- rep(0,ln_x)
    #refill_2 <- rep(0,ln_x)
    for(i in 1:ln_x){
      ##cat(" ###   I   ###" , i, "\n")
      
      if(i <= lg){
        z <- inp_dt[(1:ifelse(i==1,1,i-1)),]
        #z$promo_ind <- as.integer(z$promo_ind)
        promo_smry_1[i] <- ifelse(sum(z$promo_ind)==0,0,(sum(z$promo_ind)/min(lg,i)))
        promo_smry_2[i] <- mean(z$baseprice[z$promo_ind==1] - z$offer_price[z$promo_ind==1])
      }
      if(i > lg){
        z <- inp_dt[(i-lg):(i-1),]   # assuming lg >= 2
        #z$promo_ind <- as.integer(z$promo_ind)
        promo_smry_1[i] <- ifelse(sum(z$promo_ind)==0,0,(sum(z$promo_ind)/min(lg,i)))
        promo_smry_2[i] <- mean(z$baseprice[z$promo_ind==1] - z$offer_price[z$promo_ind==1])
      }
      
      # if(i > aa$start_vec[1]){
      
      #
      #
      oos_smry_1[i] <- ifelse(sum(z$dy_itm_loc_oos_ind)==0,0,(sum(z$dy_itm_loc_oos_ind)/min(lg,i)))
      sales_smry[i] <- sum(z$ttl_units)
      
      #Dividing the df z into 4 parts to see the rate of change of sales over weeks
      if(i <= int)
      {
        sales_rate_1[i] <- 0
        sales_rate_2[i] <- 0
        sales_rate_3[i] <- 0
      }
      else if(i <= (2*int))
      {
        s1 <- sum(z$ttl_units[1:int])
        s2 <- sum(z$ttl_units) - s1
        sales_rate_1[i] <- (s2-s1)/s1
        sales_rate_2[i] <- 0
        sales_rate_3[i] <- 0
      }
      else if(i <= (3*int))
      {
        s1 <- sum(z$ttl_units[1:int])
        s2 <- sum(z$ttl_units[1:(2*int)]) - s1
        s3 <- sum(z$ttl_units) - s1 - s2
        sales_rate_1[i] <- (s2-s1)/s1
        sales_rate_2[i] <- (s3-s2)/s2
        sales_rate_3[i] <- 0
      }
      else
      {
        s1 <- sum(z$ttl_units[1:int])
        s2 <- sum(z$ttl_units[1:(2*int)]) - s1
        s3 <- sum(z$ttl_units[1:(3*int)]) - s1 - s2
        s4 <- sum(z$ttl_units) - s1 - s2 - s3
        sales_rate_1[i] <- (s2-s1)/s1
        sales_rate_2[i] <- (s3-s2)/s2
        sales_rate_3[i] <- (s4-s3)/s3
      }
      
      y <- z$dy_itm_loc_oos_ind
      start_vec <- c()
      end_vec <- c()
      if(length(y)==1) {start_vec <- c(start_vec,1) ; end_vec <- c(end_vec,1)}
      if(length(y) > 1){
        for(j in 2:length(y)){
          if(y[j]==1 & y[j-1]==0) start_vec <- c(start_vec,j)
          if((y[j]==1 & y[j-1]==1 & j==2)) start_vec <- c(start_vec,(j-1))
          if((y[j]==0 & y[j-1]==1 & j==2)) start_vec <- c(start_vec,(j-1))
          
          if(y[j]==0 & y[j-1]==1 & j <= length(y)) end_vec <- c(end_vec,(j-1))
          if(y[j]==1 & y[j-1]==0 & j == length(y)) end_vec <- c(end_vec,j)
          if(y[j]==1 & y[j-1]==1 & j == length(y)) end_vec <- c(end_vec,j)
          
        }
      }
      oos_smry_2[i] <-  ifelse(length(end_vec) ==0 & length(start_vec)==0,0,max(end_vec - start_vec)+1)  
    }
    return(data.frame(oos_smry_1=oos_smry_1,oos_smry_2=oos_smry_2,promo_smry_1=promo_smry_1,
                      promo_smry_2=promo_smry_2, sales_smry = sales_smry, sales_rate_1 = sales_rate_1,
                      sales_rate_2 = sales_rate_2, sales_rate_3 = sales_rate_3))
  }
  
  #Create data for model with censoring
  create_data_for_model_new <- function(dataframe, lg)
  {
    aa <- inv_up_locate(dataframe$sum_boh_q,dataframe$dy_itm_loc_oos_ind)
    #aa$start_0 <- lag(aa$end_vec,1)
    #aa$start_0[1] <- 0
    bb <- promo_extractor(dataframe,lg)
    bb<-cbind(bb, dataframe$promo_ind)
    colbb <- dim(bb)[2]
    colnames(bb)[colbb] <- "promo_ind"
    bb<-cbind(bb, dataframe$sum_boh_q)
    colbb <- dim(bb)[2]
    colnames(bb)[colbb] <- "sum_boh_q"
    
    #abb <- cbind(bb,lag(dataframe$dy_itm_loc_oos_ind,1))
    model_ready_data <- bb[aa$start_vec,]
    model_ready_data$status <- dataframe$dy_itm_loc_oos_ind[aa$end_vec]
    model_ready_data$gap <- aa$end_vec - aa$start_vec + 1
    # new<-cbind(tail(abb,1), (dim(abb)[1] - tail(aa$end_vec,1)))
    # colnames(new)[(colbb+2)] <- "gap"
    # model_ready_data <- rbind(model_ready_data, new)
    model_ready_data <- model_ready_data[-1,]
    model_ready_data$sku <- rep(dataframe$sku[1], dim(model_ready_data)[1])
    
    return(model_ready_data)
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
  
  dist_to_weight <- function(distance)
  {
    weight <- 1/(distance+1)
    weight <- weight/sum(weight)
    return(weight)
  }
  
  appl_weight <- function(differ)
  {
    ifelse((differ > 0), 1.2, 0.8)
  }
  
  my_metric <- function(data, lev = NULL, model = NULL) 
  {
    diff <- data$pred - data$obs
    met <- sum(((abs(diff)/data$obs) * sapply(diff,appl_weight)))
    c(MYM = met)
  }
  
  wtd_metric <- function(expected, observed)
  {
    diff <- expected - observed
    cat("Exp",expected, "\n")
    cat("Obs",observed, "\n")
    cat("NA in expected ->", anyNA(expected), "\n")
    cat("NA in observed ->", anyNA(observed), "\n")
    cat("0 in expected ->", sum(expected == 0), "\n")
    met <- mean(((abs(diff)/expected) * sapply(diff,appl_weight)))
    return(met)
  }
  
  calc_metric <- function(df)
  {
    df <- as.data.frame(df)
    metric <- c()
    obs <- df[,1]
    for(i in 2:dim(df)[2])
    {
      met <- wtd_metric(expected = df[,i], observed = obs)
      metric <- c(metric,met)
      names(metric)[(i-1)] <- colnames(df)[i]
    }
    cat("calc_metric", metric, "\n")
    return(metric)
  }
  
  cat("START", "\n")
  
  # start_vec <- c()
  # end_vec <- c()
  # inp_x <- x[,"dy_itm_loc_oos_ind"]
  # inp_x <- as.data.frame(inp_x)
  # for(i in 2:length(inp_x)){
  #   if(inp_x[i]==1 & inp_x[i-1]==0) start_vec <- c(start_vec,i)
  #   if((inp_x[i-1]==1 & i==2)) start_vec <- c(start_vec,i-1)
  #   
  #   if(inp_x[i]==0 & inp_x[i-1]==1 & i <= length(inp_x)) end_vec <- c(end_vec,(i-1))
  #   if(inp_x[i]==1 & i == length(inp_x)) end_vec <- c(end_vec,i)
  # }
  # 
  # model_ready_data <- data.frame(start_vec=start_vec,end_vec=end_vec)
  
  #print(dim(model_ready_data))
  #cat(class(x), "\n")
  #x <- SparkR:::collect(x)
  #cat("END", "\n")
  
  model_ready_data <- try(create_data_for_model_new(x, 30));
  if(class(model_ready_data) == "try-error") model_ready_data <- NULL
  #if(dim(model_ready_data)[1] == 0) model_ready_data <- NULL
  if(is.null(model_ready_data))
  {
    metos <- rep(100, 4)
    names(metos) <- c("CoxPH_0.4_w","CoxPH_0.3_w","CoxPH_0.4","CoxPH_0.3")
  }
  model_ready_data <- na.omit(model_ready_data)
  #Greater than 10 because we have 10 covariates
  if(sum(model_ready_data$status) > 10)
  {
    cat("NA's ommited", "\n")
    
    train_dim <- dim(model_ready_data)[2] - 3
    
    #Robust Mahalanobis Distance for weighting
    #Doing this over 1:(train_dim - 1) because can ignore promo_ind for weighting
    weight_train <- try(dist_to_weight(Moutlier(model_ready_data[,2:(train_dim-2)]
                                                    ,plot = FALSE)$rd))

    if(class(weight_train) == "try-error") weight_train <- NULL
    cat("Assigned weights", "\n")
    
    thresh <- c(0.4,0.3)
    l <- list("lmodels" = c("CoxPH","CoxPH"), 
              "thresholds" = rep(thresh,1))
    looval_pred <- try(create_oosample_predictions_df_w(formula = Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+
                                                  promo_smry_2+sales_smry+sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q, models = l,
                                                data = model_ready_data, weight = weight_train,
                                                metric = "MYM", maximize = FALSE,
                                                trControl = trainControl(summaryFunction = my_metric)))
    cat(x$item_loc[1], "\n")
    #if(class(looval_pred) == "try-error") metos <- rep(150, 2)
    
    cat('Weighted models trained', "\n")
    
    looval_pred2 <- try(create_oosample_predictions_df(formula = Surv(gap, status) ~ promo_ind+oos_smry_1+
                                                         oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                                                       +sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q,
                                                       models = l,data = model_ready_data,
                                                       metric = "MYM", maximize = FALSE, 
                                                       trControl = trainControl(summaryFunction = my_metric)))
    
    cat('Unweighted models trained', "\n")
    cat('Combining results', "\n")
    if(class(looval_pred2) == "try-error")
    {
      metos <- rep(150, 4)
      names(metos) <- c("CoxPH_0.4_w","CoxPH_0.3_w","CoxPH_0.4","CoxPH_0.3")
    }
    else
    {
      looval_pred <- cbind(looval_pred,looval_pred2[-1])
      metos <- calc_metric(df = looval_pred)
    }
    cat('metos ready', "\n")
    #res <- data.frame(key, dim(model_ready_data)[1], stringsAsFactors = FALSE)
    #colnames(res) <- c("key", "numrows")
  }
  else if(sum(model_ready_data$status) > 1)
  {
    cat("Else 2 case", x$item_loc[1], "\n")
    looval <- data.frame(observed = model_ready_data$gap[model_ready_data$status == 1])
    looval$CoxPH_0.4_w <- rep(mean(looval$observed),dim(looval)[1])
    looval$CoxPH_0.3_w <- rep(mean(looval$observed),dim(looval)[1])
    looval$CoxPH_0.4 <- rep(mean(looval$observed),dim(looval)[1])
    looval$CoxPH_0.3 <- rep(mean(looval$observed),dim(looval)[1])
    
    metos <- calc_metric(df = looval)
  }
  else
  {
    metos <- rep(100, 4)
    cat("0 oos case", x$item_loc[1], "\n")
    names(metos) <- c("CoxPH_0.4_w","CoxPH_0.3_w","CoxPH_0.4","CoxPH_0.3")
  }
  #z <- data.frame(itm_loc=key,data_row=dim(y)[1], stringsAsFactors = FALSE)
  #return(y)
  res <- data.frame(t(c(key, metos)), stringsAsFactors = FALSE)
  colnames(res)[1] <- "Key"
  res
}

result2 <- gapplyCollect(all_itm_oos_agg, cols = all_itm_oos_agg$item_loc, newsparkfunc2)

newsparkfunc3 <- function(key,x)
{
  library(data.table)
  library(MTS, lib.loc = "/home_dir/z004189/Library/R")
  library(sqldf)
  library(dplyr)
  
  y <- data.frame(key, dim(x)[1], stringsAsFactors = FALSE)
  colnames(y) <- c("key", "numrows")
  #z <- data.frame(itm_loc=key,data_row=dim(y)[1], stringsAsFactors = FALSE)
  y
}


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------


#All different cases that can happen
# Case 1 -> All metrics are 100 means there is some problem in creation of model_ready_data.
# Case 2 -> All metrics are 150 means even after having more than 10 data points the model
#           is not able to fit.
# Case 3 -> All the 4 models give the same result that means data points are less than 10 so
#           average was taken.
# Case 4 -> Both unweighted and weighted models give the same result but different results
#           for different thresholds means that there was problem while assigning weights
#           that is, a singular matrix was encountered.