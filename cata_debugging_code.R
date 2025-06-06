library(dplyr)
library(ggplot2)
library(cowplot)

folds_stacked_predictions <- 
  readRDS("./Results/5foldCV_MetabricStackedPredictions.rds")

# Extract predictions for each method
predCoxPH <- folds_stacked_predictions %>%
  dplyr::select(starts_with("CoxPH"))
predCoxTime <- folds_stacked_predictions %>%
  dplyr::select(starts_with("CoxTime"))
predDeepHit <- folds_stacked_predictions %>%
  dplyr::select(starts_with("DeepHit"))
predDeepSurv <- folds_stacked_predictions %>%
  dplyr::select(starts_with("DeepSurv"))
predRSF <- folds_stacked_predictions %>%
  dplyr::select(starts_with("RSF"))

# Ensure that the predictions are matrices
predCoxPH <- as.matrix(predCoxPH)
predCoxTime <- as.matrix(predCoxTime)
predDeepHit <- as.matrix(predDeepHit)
predDeepSurv <- as.matrix(predDeepSurv)
predRSF <- as.matrix(predRSF)

# Calculate expected mortality for each method
expectedMort <- function(myPredictions) {
  -sum(log(as.numeric(myPredictions)))
}

emCoxPH <- apply(predCoxPH, 1, expectedMort)
emCoxTime <- apply(predCoxTime, 1, expectedMort)
emDeepHit <- apply(predDeepHit, 1, expectedMort)
emDeepSurv <- apply(predDeepSurv, 1, expectedMort)
emRSF <- apply(predRSF, 1, expectedMort)

em <- data.frame(
  CoxPH = emCoxPH,
  CoxTime = emCoxTime,
  DeepHit = emDeepHit,
  DeepSurv = emDeepSurv,
  RSF = emRSF
)

# PLot the expected mortality for each method
library(GGally)
em %>% ggpairs()

em[grepl("Fold 1", rownames(em)),] %>% ggpairs()
em[grepl("Fold 2", rownames(em)),] %>% ggpairs()
em[grepl("Fold 3", rownames(em)),] %>% ggpairs()
em[grepl("Fold 4", rownames(em)),] %>% ggpairs()
em[grepl("Fold 5", rownames(em)),] %>% ggpairs()

em$fold <- substr(rownames(em), 1, 6)

p1 <- em %>% 
  ggplot(aes(x = fold, y = CoxPH)) +
  geom_violin() +
  theme_minimal()

p2 <- em %>% 
  ggplot(aes(x = fold, y = CoxTime)) +
  geom_violin() +
  theme_minimal()

p3 <- em %>% 
  ggplot(aes(x = fold, y = DeepHit)) +
  geom_violin() +
  theme_minimal()

p4 <- em %>% 
  ggplot(aes(x = fold, y = DeepSurv)) +
  geom_violin() +
  theme_minimal()

p5 <- em %>% 
  ggplot(aes(x = fold, y = RSF)) +
  geom_violin() +
  theme_minimal()

p1 + p2 + p3 + p4 + p5 +
  plot_layout(ncol = 2)


mytimes <- as.numeric(gsub("CoxPH.", "", names(predCoxPH)))

mypredplot <- function(i, mytimes, predCoxPH, predCoxTime, 
                       predRSF, predDeepHit, predDeepSurv) {
  
  par(mfrow = c(1,2))
  
  plot(mytimes, as.vector(predCoxPH[i,]), type = "l")
  lines(mytimes, as.vector(predCoxTime[i,]), col = "red")
  lines(mytimes, as.vector(predRSF[i,]), col = "blue")
  lines(mytimes, as.vector(predDeepHit[i,]), col = "green")
  lines(mytimes, as.vector(predDeepSurv[i,]), col = "purple")
  title(paste("Survival curves for sample", i))
  legend("topright", legend = c("CoxPH", "CoxTime", "RSF", "DeepHit", "DeepSurv"),
         col = c("black", "red", "blue", "green", "purple"), lty = 1)
  
  plot(mytimes, -log(predCoxPH[i,]), type = "l")
  lines(mytimes, -log(predCoxTime[i,]), col = "red")
  lines(mytimes, -log(predRSF[i,]), col = "red")
  lines(mytimes, -log(predRSF[i,]), col = "blue")
  lines(mytimes, -log(predDeepHit[i,]), col = "green")
  lines(mytimes, -log(predDeepSurv[i,]), col = "purple")
  title(paste("-log(S(t)) for sample", i))
  
  cat("Expected mortality CoxPH", sum(-log(predCoxPH[i,])), ":\n")
  cat("Expected mortality CoxTime", sum(-log(predCoxTime[i,])), ":\n")
  cat("Expected mortality RSF", sum(-log(predRSF[i,])), ":\n")
  cat("Expected mortality DeepHit", sum(-log(predDeepHit[i,])), ":\n")
  cat("Expected mortality DeepSurv", sum(-log(predDeepSurv[i,])), ":\n")

}

mypredplot(i = 1, mytimes, 
           predCoxPH, predCoxTime, 
           predRSF, predDeepHit, predDeepSurv)

mypredplot(i = 2, mytimes, 
           predCoxPH, predCoxTime, 
           predRSF, predDeepHit, predDeepSurv)

mypredplot(i = 3, mytimes, 
           predCoxPH, predCoxTime, 
           predRSF, predDeepHit, predDeepSurv)

mypredplot(i = 4, mytimes, 
           predCoxPH, predCoxTime, 
           predRSF, predDeepHit, predDeepSurv)

plot(mytimes, -log(predCoxPH[2,]), type = "l")
lines(mytimes, -log(predCoxTime[2,]), col = "red")
lines(mytimes, -log(predRSF[2,]), col = "red")
lines(mytimes, -log(predRSF[2,]), col = "blue")
lines(mytimes, -log(predDeepHit[2,]), col = "green")
lines(mytimes, -log(predDeepSurv[2,]), col = "purple")

plot(mytimes, as.vector(predCoxPH[1,]), type = "l")
lines(mytimes, as.vector(predCoxPH[2,]), col = "red")

# Calculate expected mortality
-sum(log(as.numeric(predCoxPH[1,])))
-sum(log(as.numeric(predCoxPH[2,])))

  
  
  