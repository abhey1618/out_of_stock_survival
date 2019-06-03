library(data.table)
library(MTS)
library(sqldf)
library(dplyr)
library(survival)
library(ggplot2)
library(ggfortify)
library(randomForestSRC)


oos_raw <- fread("/Users/z004189/Documents/Abhey/survival/oos_sample_data_1.csv",header=T)
sum(oos_raw$temp_promo_weekly_sales_and_events_by_store.sku=='NULL')

str(oos_raw)
extr_col_nms <- unlist(strsplit(colnames(oos_raw),split='[.]'))
nw_col_nms <- extr_col_nms[seq(2,length(extr_col_nms),by=2)]
colnames(oos_raw) <- nw_col_nms

oos_raw$offer_price_orig <- as.numeric(oos_raw$offer_price_orig)
oos_raw$exit_price <- as.numeric(oos_raw$exit_price)

oos_raw$ind_elig_istk_out_of_stk_f <- as.numeric(oos_raw$sum_elig_istk_out_of_stk_f > 0)

#create oos flag from boh
oos_raw$dy_itm_loc_oos_ind <- 1 - as.numeric(oos_raw$sum_boh_q > 0)


oos_sorted <- oos_raw[,c("sku","co_loc_i","p_sls_d","p_week_start_date","co_loc_n","istk_appl_f","sum_elig_istk_out_of_stk_f",
                         "sum_boh_q","sum_eoh_q","circ_tpc_stack","baseprice","offer_price_orig","exit_price","circ_f","tpc_f",
                         "offer_price","nrgtn","ttl_units","ttl_units_gross", "dy_itm_loc_oos_ind")]
all_itm_oos_agg <-  sqldf("select sku,p_sls_d,p_week_start_date,co_loc_i,istk_appl_f,sum_elig_istk_out_of_stk_f,
circ_f,tpc_f,baseprice,offer_price_orig,exit_price,offer_price,nrgtn,ttl_units,ttl_units_gross,dy_itm_loc_oos_ind,
sum_boh_q,sum_eoh_q
                          from oos_sorted where istk_appl_f=1 group by sku,p_week_start_date,co_loc_i,p_sls_d")

sku_set=unique(oos_sorted$sku)

sku_set <- sku_set[sku_set!=949766]

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
                    sales_rate_2 = sales_rate_2, sales_rate_3 = sales_rate_3))
}

#Creates the data on which we can implement models
create_data_for_model <- function(dataframe, lg)
{
  aa <- run_locate(dataframe$dy_itm_loc_oos_ind)
  aa$start_0 <- lag(aa$end_vec,1)
  aa$start_0[1] <- 0
  bb <- promo_extractor(dataframe,lg)
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
  #model_ready_data$status <- c(rep(1,(length(model_ready_data$gap)-1)),0)
  
  return(model_ready_data)
}

all_itm_oos_agg$promo_ind <- as.numeric(all_itm_oos_agg$baseprice > all_itm_oos_agg$offer_price)

#Frequency of out of stocks for each item, location pair
oos_agg_all <- sqldf("select sku,co_loc_i,count(*) as freq from all_itm_oos_agg 
                     where dy_itm_loc_oos_ind = 1 group by sku,co_loc_i")

#Looking at a single item, location pair
#Item-> 949663, Location ->100
oos_agg_100_949663 <- all_itm_oos_agg[(all_itm_oos_agg$sku == 949663) & (all_itm_oos_agg$co_loc_i == 100),]
sum(oos_agg_100_949663$dy_itm_loc_oos_ind == 1)

#Verifying previous day's ending hand is next day's starting hand
sub<-oos_agg_100_949663[,c("sku","co_loc_i","p_sls_d","sum_boh_q","sum_eoh_q")]
#Surprisingly, its not at a lot of places

#Creating data on which we can use survival anlysis
model_ready_data <- create_data_for_model(oos_agg_100_949663, 30)

#CoxPH model
fit_oos_surv_model <- coxph(Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+
                              promo_smry_2, data=model_ready_data)

#OOS curve
nfit <- survfit(fit_oos_surv_model, newdata = model_ready_data)
plot(nfit$time,(1-nfit$surv[,2]), 'l', xlab = 'Time in days', 
     ylab = 'Out of stock probability', main = 'Out of stock probability curves', ylim = c(0,1))
for(i in 3:10)
{
  lines(nfit$time,(1-nfit$surv[,i]))
}

##Trying out Random survival forests
rfs <- rfsrc(Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2,
             data=model_ready_data, ntree = 100, samptype = "swr",
             seed = 18, na.action = "na.impute", nodesize = 5)
fore<-rfs$survival
funoob <- rfs$survival.oob

plot(rfs$time.interest,(1-rfs$survival.oob[1,]), 'l', xlab = 'Time in days', 
     ylab = 'Out of stock probability', main = 'Out of stock probability curves', ylim = c(0,1))
for(i in 2:10)
{
  lines(rfs$time.interest,(1-rfs$survival.oob[i,]))
}
newpred <- predict(rfs, newdata = model_ready_data[1:10,])
newpred$survival
plot(rfs$time.interest,(1-newpred$survival[1,]), 'l', xlab = 'Time in days', 
     ylab = 'Out of stock probability', main = 'Out of stock probability curves', ylim = c(0,1))
for(i in 2:9)
{
  lines(rfs$time.interest,(1-newpred$survival[i,]))
}

#The below commented part takes a lot of time, 
#The only use of it is to find item, location pair with most data.

# max = 0
# for(i in 1:dim(oos_agg_all)[1])
# {
#   if(oos_agg_all$freq[i] <= 2*max) ##In the best case for more data there will be oos occuring 
# every alternate day
#   {
#     next
#   }
#   dt <- all_itm_oos_agg[(all_itm_oos_agg$sku == oos_agg_all$sku[i]) &
#                           (all_itm_oos_agg$co_loc_i == oos_agg_all$co_loc_i[i]),]
#   aa <- run_locate(dt$dy_itm_loc_oos_ind)
#   if(dim(aa)[1] > max)
#   {
#     max <- dim(aa)[1]
#     ritem <- oos_agg_all$sku[i]
#     rloc <- oos_agg_all$co_loc_i[i]
#   }
# }

##Trying the same with the most data item, location pair
#ritem
#rloc
##We get that item -> 949663, location -> 3056

oos_agg_3056_949663 <- all_itm_oos_agg[(all_itm_oos_agg$sku == 949663) & (all_itm_oos_agg$co_loc_i == 3056),]
sum(oos_agg_3056_949663$dy_itm_loc_oos_ind == 1)

#Creating data on which we can use survival anlysis
model_ready_data <- create_data_for_model(oos_agg_3056_949663, 30)

#Random survival forests again
rfs <- rfsrc(Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2,
             data=model_ready_data, ntree = 100, samptype = "swr",
             seed = 18, na.action = "na.impute", nodesize = 5)
#Inbag
plot(rfs$time.interest,(1-rfs$survival[1,]), 'l', xlab = 'Time in days', 
     ylab = 'Out of stock probability', main = 'Out of stock probability curves', ylim = c(0,1))
for(i in 2:dim(rfs$survival)[1])
{
  lines(rfs$time.interest,(1-rfs$survival[i,]))
}

#Out of bag
plot(rfs$time.interest,(1-rfs$survival.oob[1,]), 'l', xlab = 'Time in days', 
     ylab = 'Out of stock probability', main = 'Out of stock probability curves', ylim = c(0,1))
for(i in 2:dim(rfs$survival.oob)[1])
{
  lines(rfs$time.interest,(1-rfs$survival.oob[i,]))
}
