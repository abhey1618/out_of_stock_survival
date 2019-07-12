# Out of stock analysis Main File
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
library(caret)
library(chemometrics)
library(heplots)

source("oos_code.R")
source("validation.R")
source("caretRegModels.R")

# This file can be found in "edge.bigred.target.com:/home_dir/z013v37/oos_sample_data_1.csv"
oos_raw <- fread("/Users/z004189/Documents/Abhey/survival/oos_sample_data_1.csv",header=T)

#Preprocessing

#Changing column names
extr_col_nms <- unlist(strsplit(colnames(oos_raw),split='[.]'))
nw_col_nms <- extr_col_nms[seq(2,length(extr_col_nms),by=2)]
colnames(oos_raw) <- nw_col_nms

#Changing data types
oos_raw$offer_price_orig <- as.numeric(oos_raw$offer_price_orig)
oos_raw$exit_price <- as.numeric(oos_raw$exit_price)
oos_raw$ind_elig_istk_out_of_stk_f <- as.numeric(oos_raw$sum_elig_istk_out_of_stk_f > 0)

oos_raw$dy_itm_loc_oos_ind <- 1 - as.numeric(oos_raw$sum_boh_q > 0)
oos_raw$sum_eoh_q <- oos_raw$sum_eoh_q/oos_raw$ttl_units
oos_raw$sum_boh_q <- oos_raw$sum_boh_q/oos_raw$ttl_units

#Subsetting the columns
oos_sorted <- oos_raw[,c("sku","co_loc_i","p_sls_d","p_week_start_date","co_loc_n","istk_appl_f","sum_elig_istk_out_of_stk_f",
                         "sum_boh_q","sum_eoh_q","circ_tpc_stack","baseprice","offer_price_orig","exit_price","circ_f","tpc_f",
                         "offer_price","nrgtn","ttl_units","ttl_units_gross", "dy_itm_loc_oos_ind")]

#Selecting only the reliable data
all_itm_oos_agg <-  sqldf("select sku,p_sls_d,p_week_start_date,co_loc_i,istk_appl_f,sum_elig_istk_out_of_stk_f,
circ_f,tpc_f,baseprice,offer_price_orig,exit_price,offer_price,nrgtn,ttl_units,ttl_units_gross,dy_itm_loc_oos_ind,
sum_boh_q,sum_eoh_q
                          from oos_sorted where istk_appl_f=1 group by sku,p_week_start_date,co_loc_i,p_sls_d")

all_itm_oos_agg$promo_ind <- as.numeric(all_itm_oos_agg$baseprice > all_itm_oos_agg$offer_price)
all_itm_oos_agg$item_loc <- paste(all_itm_oos_agg$sku, all_itm_oos_agg$co_loc_i, sep = '_')

#This file after preprocessing can be found at edge.bigred.target.com:/home_dir/z004189/oos_sample_data_2.csv

#Tring out the models for a item-location pair

#The below commented part takes a lot of time, 
#The only use of it is to find item, location pair with most data.

# oos_agg_all <- sqldf("select sku,co_loc_i,count(*) as freq from all_itm_oos_agg 
#                      where dy_itm_loc_oos_ind = 1 group by sku,co_loc_i")
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


#MODELLING
#Location -> 3056, Item -> 949663
oos_agg_3056_949663 <- all_itm_oos_agg[(all_itm_oos_agg$sku == 949663) & (all_itm_oos_agg$co_loc_i == 3056),]
sum(oos_agg_3056_949663$dy_itm_loc_oos_ind == 1)

#Creating data on which we can use survival anlysis
model_ready_data <- create_data_for_model(oos_agg_3056_949663, 30)
model_ready_data <- na.omit(model_ready_data)

#CoxPH model
fit_oos_surv_model <- coxph(Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+
                              promo_smry_2+sales_smry+sales_rate_1+sales_rate_2+sales_rate_3,
                            data=model_ready_data)

#OOS curve
nfit <- survfit(fit_oos_surv_model, newdata = model_ready_data)
plot(nfit$time,(1-nfit$surv[,1]), 'l', xlab = 'Time in days', 
     ylab = 'Out of stock probability', main = 'Out of stock probability curves', ylim = c(0,1))
for(i in 2:dim(model_ready_data)[1])
{
  if(model_ready_data$status[i] == 1) lines(nfit$time,(1-nfit$surv[,i]))
}


#Random survival forests again
rfs <- rfsrc(Surv(gap, status) ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+
               sales_smry+sales_rate_1+sales_rate_2+sales_rate_3, 
             data=model_ready_data, ntree = 100, samptype = "swr", importance = TRUE,
             seed = 18, na.action = "na.impute", nodesize = 3)
#Inbag
plot(rfs$time.interest,(1-rfs$survival[1,]), 'l', xlab = 'Time in days', 
     ylab = 'Out of stock probability', main = 'Out of stock probability curves', ylim = c(0,1), xlim = c(0,50))
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

#Trying other caret models and giving weights for training

#Gives us the number of covariates
train_dim <- dim(model_ready_data)[2] - 3

#Robust Mahalanobis Distance for weighting
#Doing this over 1:(train_dim - 1) because can ignore promo_ind for weighting
weight_train <- dist_to_weight(Moutlier(model_ready_data[,3:(train_dim-1)]
                                        ,plot = FALSE)$rd)

newm <- train(x = model_ready_data[,1:train_dim], y = model_ready_data$gap, method = "bagEarth", 
              weights = weight_train, metric = "MYM", maximize = FALSE, 
              trControl = trainControl(summaryFunction = my_metric))
exp <- predict(newm, newdata = model_ready_data[,1:train_dim])
exp <- ceiling(exp)
View(cbind(exp, model_ready_data$gap))

newm2 <- train(x = model_ready_data[,1:train_dim], y = model_ready_data$gap, method = "rf",
               metric = "RMSE", maximize = FALSE)
exp2 <- predict(newm2, newdata = model_ready_data[,1:train_dim])
exp2 <- ceiling(exp2)
View(cbind(exp2, model_ready_data$gap))

##Validation
cpr1 <- leave_one_out(data = model_ready_data, model = "CoxPH", formula = Surv(gap, status) ~ 
                        promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                      +sales_rate_1+sales_rate_2+sales_rate_3)
exp1 <- oospredict(prob_matrix = cpr1, threshold = 0.7)
cpr2 <- leave_one_out(data = model_ready_data, model = "RandomSurvivalForest", formula = Surv(gap, status) ~ 
                        promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                      +sales_rate_1+sales_rate_2+sales_rate_3,weights = weight_train)
exp2 <- oospredict(prob_matrix = cpr2, threshold = 0.7)
View(cbind(exp1,exp2,model_ready_data$gap))

##Censoring
model_ready_data <- create_data_for_model_new(oos_agg_3056_949663, 30)
model_ready_data <- na.omit(model_ready_data)

pr1 <- leave_one_out_cnsr(data = model_ready_data, model = "CoxPH", formula = Surv(gap, status) ~ 
                            promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                          +sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q, plot = TRUE)
exp1 <- oospredict(prob_matrix = pr1, threshold = 0.7)
pr2 <- leave_one_out_cnsr(data = model_ready_data, model = "RandomSurvivalForest", formula = Surv(gap, status) ~ 
                            promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                          +sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q)#,weights = weight_train)
exp2 <- oospredict(prob_matrix = pr2, threshold = 0.7)
View(cbind(exp1,exp2,model_ready_data[model_ready_data$status == 1,]$gap))

##Creating dataframe for results with censoring
thresh <- c(0.4,0.3) #Found through cross-validation
l <- list("lmodels" = c("RandomSurvivalForest","RandomSurvivalForest",
                        "CoxPH","CoxPH"), 
          "thresholds" = rep(thresh,2))

weight_train <- dist_to_weight(Moutlier(model_ready_data[,1:(train_dim-3)]
                                        ,plot = FALSE)$rd)

looval_pred <- create_oosample_predictions_df_w(formula = Surv(gap, status) ~ promo_ind+oos_smry_1+
                                                  oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                                                +sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q,
                                                formula2 = gap ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+
                                                  sales_smry+sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q,
                                                models = l,data = model_ready_data, weights = weight_train,
                                                metric = "MYM", maximize = FALSE, 
                                                trControl = trainControl(summaryFunction = my_metric))


looval_pred2 <- create_oosample_predictions_df(formula = Surv(gap, status) ~ promo_ind+oos_smry_1+
                                                 oos_smry_2+promo_smry_1+promo_smry_2+sales_smry
                                               +sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q,
                                               formula2 = gap ~ promo_ind+oos_smry_1+oos_smry_2+promo_smry_1+promo_smry_2+
                                                 sales_smry+sales_rate_1+sales_rate_2+sales_rate_3+sum_boh_q,
                                               models = l,data = model_ready_data,
                                               metric = "MYM", maximize = FALSE, 
                                               trControl = trainControl(summaryFunction = my_metric))

looval_pred <- cbind(looval_pred,looval_pred2[-1])