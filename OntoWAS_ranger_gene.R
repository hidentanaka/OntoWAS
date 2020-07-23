library(readr)
library(dplyr)
library(assertr)
library(rsample)

require(data.table)
require(ranger)



data.dt <- fread("inputfile_gene.csv")
head(data.dt)


#K-fold = 10
train_folds <- vfold_cv(data.dt, v = 5)
train_folds

train1.dt <- analysis(train_folds$splits[[1]])
test1.dt <- assessment(train_folds$splits[[1]])

train2.dt <- analysis(train_folds$splits[[2]])
test2.dt <- assessment(train_folds$splits[[2]])

train3.dt <- analysis(train_folds$splits[[3]])
test3.dt <- assessment(train_folds$splits[[3]])

train4.dt <- analysis(train_folds$splits[[4]])
test4.dt <- assessment(train_folds$splits[[4]])

train5.dt <- analysis(train_folds$splits[[5]])
test5.dt <- assessment(train_folds$splits[[5]])


#data = 1#

ranger_model <- ranger(phenotype ~ ., data = train1.dt, num.trees = 1000, importance = 'permutation')


PIMP <- importance_pvalues(ranger_model, method = "altmann", formula = phenotype ~ ., data = train1.dt)
#iteration=5000
#importance_pvalues(ranger_model, method = "altmann",num.permutations = 5000, formula = phenotype ~ ., data = train.dt)

#predict
ranger_pred <- predict(ranger_model, data=test1.dt)
#predict_list
ranger_pred$predictions
#R^2 score
actual <- test1.dt$phenotype
predicted <- ranger_pred$predictions
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
print(R2)


sink('Phenotypename_gene_test1.txt', append = TRUE)
print (PIMP)
sink()

#data = 2#

ranger_model <- ranger(phenotype ~ ., data = train1.dt, num.trees = 1000, importance = 'permutation')


PIMP <- importance_pvalues(ranger_model, method = "altmann", formula = phenotype ~ ., data = train1.dt)
#iteration=5000
#importance_pvalues(ranger_model, method = "altmann",num.permutations = 5000, formula = phenotype ~ ., data = train.dt)

#predict
ranger_pred <- predict(ranger_model, data=test1.dt)
#predict_list
ranger_pred$predictions
#R^2 score
actual <- test1.dt$phenotype
predicted <- ranger_pred$predictions
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
print(R2)


sink('Phenotypename_gene_test2.txt', append = TRUE)
print (PIMP)
sink()

#data = 3#

ranger_model <- ranger(phenotype ~ ., data = train1.dt, num.trees = 1000, importance = 'permutation')


PIMP <- importance_pvalues(ranger_model, method = "altmann", formula = phenotype ~ ., data = train1.dt)
#iteration=5000
#importance_pvalues(ranger_model, method = "altmann",num.permutations = 5000, formula = phenotype ~ ., data = train.dt)

#predict
ranger_pred <- predict(ranger_model, data=test1.dt)
#predict_list
ranger_pred$predictions
#R^2 score
actual <- test1.dt$phenotype
predicted <- ranger_pred$predictions
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
print(R2)


sink('Phenotypename_gene_test3.txt', append = TRUE)
print (PIMP)
sink()


#data = 4#

ranger_model <- ranger(phenotype ~ ., data = train1.dt, num.trees = 1000, importance = 'permutation')


PIMP <- importance_pvalues(ranger_model, method = "altmann", formula = phenotype ~ ., data = train1.dt)
#iteration=5000
#importance_pvalues(ranger_model, method = "altmann",num.permutations = 5000, formula = phenotype ~ ., data = train.dt)

#predict
ranger_pred <- predict(ranger_model, data=test1.dt)
#predict_list
ranger_pred$predictions
#R^2 score
actual <- test1.dt$phenotype
predicted <- ranger_pred$predictions
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
print(R2)


sink('Phenotypename_gene_test4.txt', append = TRUE)
print (PIMP)
sink()

#data = 5#

ranger_model <- ranger(phenotype ~ ., data = train1.dt, num.trees = 1000, importance = 'permutation')


PIMP <- importance_pvalues(ranger_model, method = "altmann", formula = phenotype ~ ., data = train1.dt)
#iteration=5000
#importance_pvalues(ranger_model, method = "altmann",num.permutations = 5000, formula = phenotype ~ ., data = train.dt)

#predict
ranger_pred <- predict(ranger_model, data=test1.dt)
#predict_list
ranger_pred$predictions
#R^2 score
actual <- test1.dt$phenotype
predicted <- ranger_pred$predictions
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))
print(R2)


sink('Phenotypename_gene_test5.txt', append = TRUE)
print (PIMP)
sink()
