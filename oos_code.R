library(data.table)
library(MTS)
library(sqldf)
library(dplyr)
library(survival)
library(ggplot2)
library(ggfortify)

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
  #lg <- 21
  ln_x = dim(inp_dt)[1]
  oos_smry_1 <- rep(0,ln_x)
  oos_smry_2 <- rep(0,ln_x)
  promo_smry_1 <- rep(0,ln_x)
  promo_smry_2 <- rep(0,ln_x)
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
    #dt <- run_locate(x$dy_itm_loc_oos_ind)
    
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
  return(data.frame(oos_smry_1=oos_smry_1,oos_smry_2=oos_smry_2,promo_smry_1=promo_smry_1,promo_smry_2=promo_smry_2))
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
aa <- run_locate(oos_agg_100_949663$dy_itm_loc_oos_ind)
bb <- promo_extractor(oos_agg_100_949663,5)
bb<-cbind(bb, oos_agg_100_949663$promo_ind)
colnames(bb)[5] <- "promo_ind"

abb <- cbind(oos_agg_100_949663$dy_itm_loc_oos_ind,bb)
model_ready_data <- abb[aa$start_vec,]

model_ready_data$gap <- aa$end_vec - aa$start_vec +1 
model_ready_data$status <- c(rep(1,(length(model_ready_data$gap)-1)),0)

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
library(randomForestSRC)
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

##Trying the same with the most data item, location pair
max(oos_agg_all[,'freq'])
oos_agg_all[oos_agg_all$freq == 279,]
##We get that item -> 1567262, location -> 1720

oos_agg_1720_1567262 <- all_itm_oos_agg[(all_itm_oos_agg$sku == 1567262) & (all_itm_oos_agg$co_loc_i == 1720),]
sum(oos_agg_1720_1567262$dy_itm_loc_oos_ind == 1)
