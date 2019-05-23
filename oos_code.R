library(data.table)
library(MTS)
library(sqldf)
library(dplyr)
library(survival)

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
                         "offer_price","nrgtn","ttl_units","ttl_units_gross")]
all_itm_oos_agg <-  sqldf("select sku,p_sls_d,p_week_start_date,co_loc_i,istk_appl_f,sum_elig_istk_out_of_stk_f,
circ_f,tpc_f,baseprice,offer_price_orig,exit_price,offer_price,nrgtn,ttl_units,ttl_units_gross, count(*) as oos_agg
                          from oos_sorted where sum_eoh_q<= 0 and  istk_appl_f=1 group by sku,p_week_start_date,co_loc_i,p_sls_d")

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


