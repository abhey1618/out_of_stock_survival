library(data.table)
library(MTS)
library(sqldf)
library(dplyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(randomForestSRC)
library(coxme)
library(quantreg)

#Function to locate runs
run_locate <- function(x){
  start_vec <- c()
  end_vec <- c()
  for(i in 2:length(x)){
    if(x[i]==1 & x[i-1]==0) start_vec <- c(start_vec,i)
    if((x[i-1]==1 & i==2)) start_vec <- c(start_vec,i-1)
    
    if(x[i]==0 & x[i-1]==1 & i <= length(x)) end_vec <- c(end_vec,(i-1))
    if(x[i]==1 & i == length(x)) end_vec <- c(end_vec,i)
  }
  return(data.frame(start_vec=start_vec,end_vec=end_vec))
}

#Function to locate censoring
inv_up_locate <- function(x,y){
  start_vec <- c()
  end_vec <- c()
  mark1 <- 0
  mark2 <- 0
  for(i in 2:length(x)){
    #mark1 is to remember start of an oos
    #mark2 is to see if an oos is still going on or not
    if(y[i] == 1 & mark2 == 0)
    {
      mark1 <- i
      mark2 <- 1
    }
    if((x[i]>0 & i==2)) start_vec <- c(start_vec,1)
    if(x[i] > x[i-1] & y[i] == 0)
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
    
    #if(x[i]==0 & x[i-1]==1 & i <= length(x)) end_vec <- c(end_vec,(i-1))
    #if(x[i]==1 & i == length(x)) end_vec <- c(end_vec,i)
  }
  end_vec <- c(end_vec,length(x))
  if(end_vec[1] == 0)
  {
    start_vec <- start_vec[-1]
    end_vec <- end_vec[-1]
  }
  if(end_vec[1] < start_vec[1])
  {
    end_vec <- end_vec[-1]
  }
  #return(c(length(start_vec),length(end_vec)))
  return(data.frame(start_vec=start_vec,end_vec=end_vec))
}

#Function to generate covariates (Without censoring case)
promo_extractor_old <- function(inp_dt,lg){
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
  refill_1 <- rep(0,ln_x)
  refill_2 <- rep(0,ln_x)
  for(i in 1:ln_x){
    #cat(" ###   I   ###" , i, "\n")
    
    if(i <= lg){
      x <- inp_dt[(1:ifelse(i==1,1,i-1)),]
      promo_smry_1[i] <- ifelse(sum(x$promo_ind)==0,0,(sum(x$promo_ind)/min(lg,i)))
      promo_smry_2[i] <- mean(x$baseprice[x$baseprice!=x$exit_price] - x$exit_price[x$baseprice!=x$exit_price])
    }
    if(i > lg){
      x <- inp_dt[(i-lg):(i-1),]   # assuming lg >= 2
      promo_smry_1[i] <- ifelse(sum(x$promo_ind)==0,0,(sum(x$promo_ind)/min(lg,i)))
      promo_smry_2[i] <- mean(x$baseprice[x$baseprice!=x$exit_price] - x$exit_price[x$baseprice!=x$exit_price])
    }
    
    # if(i > aa$start_vec[1]){
    
    #
    #
    oos_smry_1[i] <- ifelse(sum(x$dy_itm_loc_oos_ind)==0,0,(sum(x$dy_itm_loc_oos_ind)/min(lg,i)))
    sales_smry[i] <- sum(x$ttl_units)
    
    fr <- 0
    mx_rfll <- 0
    temp <- 0
    if(i >= 3)
    {
      boh <- x$sum_boh_q
      for(j in 2:length(boh))
      {
        if(boh[j] > boh[(j-1)])
        {
          fr <- fr+1
          if(j - temp > mx_rfll)
          {
            mx_rfll <- j - temp
            temp <- j
          }
        }
      }
    }
    refill_1[i] <- fr
    refill_2[i] <- mx_rfll
    #Dividing the df x into 4 parts to see the rate of change of sales over weeks
    if(i <= int)
    {
      sales_rate_1[i] <- 0
      sales_rate_2[i] <- 0
      sales_rate_3[i] <- 0
    }
    else if(i <= (2*int))
    {
      s1 <- sum(x$ttl_units[1:int])
      s2 <- sum(x$ttl_units) - s1
      sales_rate_1[i] <- (s2-s1)/s1
      sales_rate_2[i] <- 0
      sales_rate_3[i] <- 0
    }
    else if(i <= (3*int))
    {
      s1 <- sum(x$ttl_units[1:int])
      s2 <- sum(x$ttl_units[1:(2*int)]) - s1
      s3 <- sum(x$ttl_units) - s1 - s2
      sales_rate_1[i] <- (s2-s1)/s1
      sales_rate_2[i] <- (s3-s2)/s2
      sales_rate_3[i] <- 0
    }
    else
    {
      s1 <- sum(x$ttl_units[1:int])
      s2 <- sum(x$ttl_units[1:(2*int)]) - s1
      s3 <- sum(x$ttl_units[1:(3*int)]) - s1 - s2
      s4 <- sum(x$ttl_units) - s1 - s2 - s3
      sales_rate_1[i] <- (s2-s1)/s1
      sales_rate_2[i] <- (s3-s2)/s2
      sales_rate_3[i] <- (s4-s3)/s3
    }
    
    y <- x$dy_itm_loc_oos_ind
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
      #cat(" ###   i  ###",i, "  ### len  diff  : startvc - endvec ##", length(start_vec) - length(end_vec),"\n")
      #cat(" $  $  $  i  $  $  $ ", i, "###", end_vec, "###", start_vec,  "\n")
      
      #oos_smry_2[i] <-  ifelse(length(end_vec) ==0 & length(end_vec)==0,0,max(end_vec - start_vec)+1)
      #cat(" ### i  ### ",i,"stvec", start_vec," endvec" , end_vec, "\n")
      
      #}
      #cat(" ### start_vec[1]  ### ",aa$start_vec[length(start_vec)]," ###   I   ###" , i, "\n")
      #cat(" ### end_vec[1]  ### ",aa$end_vec[length(end_vec)]," ###   I   ###" , i, "\n")
      
      
      # cat("###  K ",k,"\n")
      
      # if(i > lg){
      #   x <- inp_dt[(i-lg):i,]
      #   oos_smry_1[i] <- ifelse(sum(x$oos_agg)==0,0,(sum(x$oos_agg)/i))
      #   tgt_ind <- which.min((aa$start_vec-i)[aa$start_vec-i > 0])
      #   oos_smry_2[i] <- max(aa$end_vec[1:tgt_ind] - aa$start_vec[1:tgt_ind])
    }
    oos_smry_2[i] <-  ifelse(length(end_vec) ==0 & length(start_vec)==0,0,max(end_vec - start_vec)+1)  
  }
  return(data.frame(oos_smry_1=oos_smry_1,oos_smry_2=oos_smry_2,promo_smry_1=promo_smry_1,
                    promo_smry_2=promo_smry_2, sales_smry = sales_smry, sales_rate_1 = sales_rate_1,
                    sales_rate_2 = sales_rate_2, sales_rate_3 = sales_rate_3, refill_1 = refill_1,
                    refill_2 = refill_2))
}

#Generates covariates(With censoring case)
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

#Creates the data on which we can implement models
create_data_for_model <- function(dataframe, lg)
{
  aa <- run_locate(dataframe$dy_itm_loc_oos_ind)
  aa$start_0 <- lag(aa$end_vec,1)
  aa$start_0[1] <- 0
  bb <- promo_extractor_old(dataframe,lg)
  bb<-cbind(bb, dataframe$promo_ind)
  colbb <- dim(bb)[2]
  colnames(bb)[colbb] <- "promo_ind"
  
  abb <- cbind(bb,dataframe$dy_itm_loc_oos_ind)
  colnames(abb)[(colbb+1)] <- "status"
  model_ready_data <- abb[aa$start_vec,]
  model_ready_data$gap <- aa$start_vec - aa$start_0 -1
  new<-cbind(tail(abb,1), (dim(abb)[1] - tail(aa$end_vec,1)))
  colnames(new)[(colbb+2)] <- "gap"
  model_ready_data <- rbind(model_ready_data, new)
  model_ready_data$sku <- rep(dataframe$sku[1], dim(model_ready_data)[1])
  
  return(model_ready_data)
}

#creates data based on censoring according to inv_up_locate
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
