library(caret)
library(chemometrics)
library(heplots)

dist_to_weight <- function(distance)
{
  weight <- 1/(distance+1)
  weight <- weight/sum(weight)
  return(weight)
}

appl_weight <- function(x)
{
  ifelse((x > 0), 1.2, 0.8)
}

my_metric <- function(data, lev = NULL, model = NULL) 
{
  diff <- data$pred - data$obs
  met <- sum(((abs(diff)/data$obs) * sapply(diff,appl_weight)))
  c(MYM = met)
}