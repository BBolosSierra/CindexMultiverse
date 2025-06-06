library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(cowplot)
library(GGally)

# Load predictions across all folds
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

# Add sample identifiers and fold information
em <- em %>% 
  rownames_to_column(var = "Sample") %>%
  mutate(fold = substr(Sample, 1, 6)) %>%
  arrange(fold)

# Calculate the number of infinite values per fold and method
em %>%
  group_by(fold) %>%
  summarise(
    CoxPH_inf = sum(is.infinite(CoxPH)),
    CoxTime_inf = sum(is.infinite(CoxTime)),
    DeepHit_inf = sum(is.infinite(DeepHit)),
    DeepSurv_inf = sum(is.infinite(DeepSurv)),
    RSF_inf = sum(is.infinite(RSF))
  )

## As it can be seen above, there is a high number of infinite values for the 
## CoxTime and DeepSurv methods. These correspond to cases in which the survival
## function was predicted to be exactly equal to zero
  
# Plot the expected mortality for each method
em %>% 
  dplyr::select(CoxPH, CoxTime, DeepHit, DeepSurv, RSF) %>%
  ggpairs()

# As the plot above suggested multimodality for DeepHit, we now stratify by fold
for(i in unique(em$fold)) {
  cat("Fold:", i, "\n")
  p <- em %>% filter(fold == i) %>% 
    dplyr::select(CoxPH, CoxTime, DeepHit, DeepSurv, RSF) %>%
    ggpairs() +
    ggtitle(paste("Pairwise plots for fold", i))
  print(p)
}

# Plot the distribution of expected mortality for each method across folds
em_long <- em %>%
  pivot_longer(cols = c(CoxPH, CoxTime, DeepHit, DeepSurv, RSF), 
               names_to = "Method", values_to = "ExpectedMortality")

em_long %>%
  ggplot(aes(x = fold, y = ExpectedMortality, fill = fold)) +
  geom_violin() +
  facet_wrap(~ Method) +
  theme_minimal() +
  labs(title = "Expected mortality distribution by method and fold")

# Now we will plot the survival curves for a few samples
mytimes <- as.numeric(gsub("CoxPH.", "", colnames(predCoxPH)))

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

  
  
  