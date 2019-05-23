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