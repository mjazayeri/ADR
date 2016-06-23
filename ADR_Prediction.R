setwd("/Users/Mehdi/Desktop/Bioinformatics/FinalProject")
library(impute)
library(FSelector)
library(class)
library(e1071)
library(rpart)
#----------------------------- DATA LOADING AND PREPERATION -----------------------------
gene.data <- read.csv("Gene_Expression_681x978.csv", header=T)
dimnames(gene.data)[[1]] <- gene.data[, 1]
gene.data <- gene.data[, -1]

maccs.data <- read.csv("MACCS_681x166.csv", header = T)
maccs.data <- maccs.data[, -1]
maccs.valid <- which(!apply(maccs.data, 2, FUN = function(x) {all(x==0)}))
maccs.data <- maccs.data[, maccs.valid]

adr.data <- read.csv("FAERS_ADR_681x9404.csv", header = T)
dimnames(adr.data)[[1]] <- adr.data[, 1]
adr.data <- adr.data[, -1]
dist <- sapply(1:ncol(adr.data), function(x) { sum(adr.data[, x] == 0) / nrow(adr.data) })

# remove ADRs that contains 75% of 0's 
adr.best.ind <- which(dist <  0.75)
adr.data <- adr.data[, adr.best.ind]
#----------------------------- FEATURE SELECTION FUNCTIONS -----------------------------

# selects genes with difference higher than a threshold between max and min
select_predictive_features <- function(threshold) {
  diffs <- sapply(1:ncol(gene.data), function(x) {max(gene.data[, x]) - min(gene.data[, x])})
  return(which(diffs > threshold))
}

adr <- NULL
correlation_feature_selection <- function(test.row.ind, adr_col, threshold) {
  
  adr <<- adr.data[-test.row.ind, adr_col]
  cors <- linear.correlation(adr~., data = final.data[-test.row.ind, ])
  features.ind <- which(cors$attr_importance > threshold)
  if(length(features.ind) < 20) {
    top20 <- order(cors$attr_importance, decreasing = T)[1:20]
    return (top20)
  }
  return (features.ind)
}

#test set indeces are smaller, 
#so it's better to pass them as argument rather than train set indeces
t.test_feature_selection <-function(test.row.ind, adr_col, threshold) {
  zeros <- which(adr.data[, adr_col] == 0)
  ones <- which(adr.data[, adr_col] == 1)
  
  p.vect <- apply(final.data[-test.row.ind, ], 2, function(x) { t.test(x[zeros], x[ones])$p.value })
  features.ind <- which(p.vect < threshold)
  if(length(features.ind) < 20) {
    top20 <- order(p.vect, decreasing = T)[1:20]
    return (top20)
  }
  return (features.ind)
}

# generates a formula. ex: Y ~ X3 + X4 + X7
get_formula <- function(label.name, col.indeces) {
  variables <- paste(colnames(final.data[, col.indeces]), collapse = "+")
  formula <- as.formula(paste(label.name, variables))
  
  return(formula)
}

#----------------------------- CLASSIFICATION FUNCTIONS -----------------------------
adr <- NULL
get_kth_fold <- function(sample, i){
  
  fold.size <- ceiling(nrow(final.data) / 5)
  start <- 1 + ((i-1) * fold.size)
  #fold.size = end - start + 1 -> end = fold.size + start -1
  end <- start + fold.size - 1;
  #in case end was outside the boundry
  end <- min(nrow(final.data), end)
  
  
  test.ind <- sample[start : end]
  return (test.ind)
}
evaluate_knn <- function (test, features, adr_col) {

  ret.knn <- knn(train = final.data[-test, features], 
                test = final.data[test, features], 
                cl = adr.data[-test, adr_col], 
                k = 10)
  #error rate
  return( sum(ret.knn != adr.data[test, adr_col]) / length(test) )
}
evaluate_svm <- function (test, features, adr_col) {

  adr <<- adr.data[-test, adr_col]
  formula <- get_formula("adr~", features)
  svm.model <- svm(formula, data = final.data[-test, features])
  
  ret.svm = predict(svm.model, final.data[test, features])
  ret.svm = round(ret.svm)
  
  #error rate
  return( sum(ret.svm != adr.data[test, adr_col]) / length(test) )
}
evaluate_decisiontree <- function (test, features, adr_col) {
  
  adr <<- adr.data[-test, adr_col]
  formula <- get_formula("adr~", features)
  
  #controls the growth of the tree
  ctrl <- rpart.control(minsplit=2, minbucket=1, cp=0.0001)
  fit.model <- rpart(formula = formula, data = final.data[-test, features], method = "anova", control = ctrl) 
  
  #prune the descision tree to avoid overfitting
  # tree.model <- prune(fit.model, cp = fit.model$cptable[which.min(fit.model$cptable[,"xerror"]),"CP"])
  
  # ret.tree <- predict(tree.model, final.data[test, features])
  ret.tree <- predict(fit.model, final.data[test, features])
  ret.tree <- round(ret.tree)
  #error rate
  return( sum(ret.tree != adr.data[test, adr_col]) / length(test) )
}
evaluate_forwardSelection <- function (test, adr_col) {
  y <- adr.data[-test, adr_col]
  fit.null <- glm(y ~ 1, data = final.data[-test, ], family = "binomial")
  fit.full <- glm(y ~ ., data = final.data[-test, ], family = "binomial")
  invisible(fit.best <- step(fit.null, scope = list(upper = fit.full), data = final.data[-test, ], direction = "forwar"))
  ret.forward <- predict(fit.best, final.data[test, ])
  ret.forward <- round(ret.svm)
  
  return ( sum(ret.forward != adr.data[test, adr_col]) / length(test) )
}
predict_adverse_drug_reaction <- function(rounds, adr.col.ind) {
  knn.null <- NULL; knn.cor <- NULL; knn.t.test <- NULL
  svm.null <- NULL; svm.cor <- NULL; svm.t.test <- NULL
  dtree.null <- NULL; dtree.cor <- NULL; dtree.t.test <- NULL
  forward_selection <- NULL
  
  print(paste("adr_col:", adr.col.ind))
  #repeat each 10-fold cv 5 times
  for(i in 1 : rounds) {
    #create a random sequence of all indecs
    cv.sample.ind <- sample(1:nrow(final.data))
    
    #10-fold cross validation for each of 8 combination
    for(k in 1 : 5) {
      test.ind <- get_kth_fold(cv.sample.ind, k)
      
      # column index of selected features
      # no feature selection => all possible columns
      # features.null <- 1 : ncol(final.data)
      features.cor<- correlation_feature_selection(test.ind, adr.col.ind, 0.12) 
      features.t.test <- t.test_feature_selection(test.ind, adr.col.ind, 0.05)
      
      # knn combinations
      # knn.null <- c(evaluate_knn(test.ind, features.null, adr.col.ind), knn.null)
      knn.cor <- c(knn.cor, evaluate_knn(test.ind, features.cor, adr.col.ind))
      knn.t.test <- c(knn.t.test, evaluate_knn(test.ind, features.t.test, adr.col.ind))
      
      # svm combinations 
      # svm.null <- c(evaluate_svm(test.ind, features.null, adr.col.ind), svm.null)
      svm.cor <- c(evaluate_svm(test.ind, features.cor, adr.col.ind), svm.cor)
      svm.t.test <- c(evaluate_svm(test.ind, features.t.test, adr.col.ind), svm.t.test)
      
      # decision tree combinations
      # dtree.null <- c(evaluate_decisiontree(test.ind, features.null, adr.col.ind), dtree.null)
      dtree.cor <- c(evaluate_decisiontree(test.ind, features.cor, adr.col.ind), dtree.cor)
      dtree.t.test <- c(evaluate_decisiontree(test.ind, features.t.test, adr.col.ind), dtree.t.test)
      
      forward_selection <- c(evaluate_forwardSelection(test.ind, adr.col.ind), forward_selection)
    }
  }
  # take the average of all 50 runs for each combination and create a row
  res.run <- data.frame(knn_cor = mean(knn.cor), knn_t.test  = mean(knn.t.test),
                        svm_cor = mean(svm.cor), svm_t.test = mean(svm.t.test),
                        dtree_cor = mean(dtree.cor), dtree_t.test = mean(dtree.t.test),
                        forward_selection = mean(forward.selection)
  )
  #add the row as a result of evaluations for adverse drug effect of adr.col.ind
  # res.final <<- rbind(res.final, res.run)
  write.table(res.run, file="232ADR-Predictions.txt", sep = "\t", col.names = F, row.names = F, append = T)
}

#----------------------------- MAIN SESSION -----------------------------
best.genes <- select_predictive_features(0.2)
final.data <- cbind(gene.data[, best.genes], maccs.data)

#main loop that, column by column, predicts adverse drug reaction
sapply( 1:ncol(adr.data), function(x) predict_adverse_drug_reaction(rounds = 1, adr.col.ind = x) )


#----------------------------- TEST SESSION -----------------------------

adr_col <- 2
res_svm <- function (test, features, adr_col) {
  adr <<- adr.data[-test, adr_col]
  formula <- get_formula("adr~", features)
  svm.model <- svm(formula, data = final.data[-test, features])
  
  ret.svm = predict(svm.model, final.data[test, features])
  return(ret.svm)
}
res_decisiontree <- function (test, features, adr_col) {
  
  adr <<- adr.data[-test, adr_col]
  formula <- get_formula("adr~", features)
  
  #controls the growth of the tree
  ctrl <- rpart.control(minsplit=2, minbucket=1, cp=0.0001)
  fit.model <- rpart(formula = formula, data = final.data[-test, features], method = "anova", control = ctrl) 
  
  ret.tree <- predict(fit.model, final.data[test, features])
  return(ret.tree)
}
ret.svm.cor <- NULL
ret.svm.t <- NULL

cv.sample.ind <- sample(1:nrow(final.data))
for(i in 1:5) {
  test.ind <- get_kth_fold(cv.sample.ind, i)
  features.cor<- correlation_feature_selection(test.ind, adr_col, 0.12) 
  features.t.test <- t.test_feature_selection(test.ind, adr_col, 0.05)
  
  ret.svm.t <- c(ret.svm.t, res_svm(test.ind, features.t.test, adr_col))
  ret.svm.cor <- c(ret.svm.cor, res_svm(test.ind, features.cor, adr_col))
}

rocobj1 <- plot.roc(adr.data[cv.sample.ind, adr_col], ret.svm.t, percent=TRUE, col="black")  
rocobj2 <- lines.roc(adr.data[cv.sample.ind, adr_col], ret.svm.cor, percent=TRUE, col="red")  
