library(caret)
library(chemometrics)

dist_to_weight <- function(distance)
{
  weight <- 1/distance
  weight <- weight/sum(weight)
  return(weight)
}

weight_train <- dist_to_weight(Moutlier(model_ready_data[,c(1:5,8)], plot = FALSE)$rd)

appl_weight <- function(x)
{
  ifelse((x > 0), 0.6, 0.4)
}

my_metric <- function(data, lev = NULL, model = NULL) 
{
  diff <- data$pred - data$obs
  met <- sum(((diff^2) * sapply(diff,appl_weight)))
  c(MYM = met)
}

model_ready <- model_ready_data
train_dim <- dim(model_ready)[2] - 2
newm <- train(x = model_ready[,1:train_dim], y = model_ready$gap, method = "rf", weights = weight_train,
              metric = "MYM", maximize = FALSE, trControl = trainControl(summaryFunction = my_metric))
exp <- predict(newm, newdata = model_ready[,1:train_dim])
View(cbind(exp, model_ready$gap))

newm2 <- train(x = model_ready_data[,1:train_dim], y = model_ready_data$gap, method = "earth",
               metric = "RMSE", maximize = FALSE)
exp2 <- predict(newm2, newdata = model_ready_data[,1:train_dim])
View(cbind(exp2, model_ready_data$gap))