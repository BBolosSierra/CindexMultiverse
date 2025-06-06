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

n.sel <- 100
ids1 <- em$Sample[em$fold == "Fold 1"][1:n.sel]
ids2 <- em$Sample[em$fold == "Fold 2"][1:n.sel]
ids3 <- em$Sample[em$fold == "Fold 3"][1:n.sel]
ids4 <- em$Sample[em$fold == "Fold 4"][1:n.sel]
ids5 <- em$Sample[em$fold == "Fold 5"][1:n.sel]

plot(mytimes, predDeepHit[ids1[1],], type = "l")
for(i in 2:n.sel) {
  lines(mytimes, predDeepHit[ids1[i],])
}
for(i in 1:n.sel) {
  lines(mytimes, predDeepHit[ids2[i],], col = "red")
#  lines(mytimes, predDeepHit[ids3[i],], col = "blue")
#  lines(mytimes, predDeepHit[ids4[i],], col = "green")
#  lines(mytimes, predDeepHit[ids5[i],], col = "purple")
}
legend("topright", legend = c("Fold 1", "Fold 2", "Fold 3", "Fold 4", "Fold 5"),
       col = c("black", "red", "blue", "green", "purple"), lty = 1)


# Maximum times per fold
folds_stacked_predictions %>%
  group_by(cv_fold) %>%
  summarise(max_time = max(test_time),
            max_obs_time = max(test_time[test_status == 1])) %>%
  arrange(cv_fold)


folds_stacked_predictions %>%
  ggplot(aes(x = as.factor(cv_fold), y = test_time, fill = as.factor(cv_fold))) +
  geom_violin() +
  theme_minimal() +
  labs(title = "Distribution of survival times by fold")

library(survival)
library(survminer)

fit <- survfit(Surv(test_time, test_status) ~ as.factor(cv_fold), 
               data = folds_stacked_predictions)
ggsurvplot(fit, data = folds_stacked_predictions, pval = FALSE, conf.int = TRUE)


# Example from Sonabend et al paper
library(survival)
library(ranger)
library(distr6)
library(mboost)
library(mlr3proba)

set.seed(20220109)

## Get data and variables
data = survival::rats
data$sex = as.integer(data$sex == "f")
train = sample(nrow(data), nrow(data) * 2/3)
test = setdiff(seq(nrow(data)), train)
test_unique_times = data$time[test]
target = Surv(test_unique_times, data$status[test])
target_train = Surv(data$time[train], data$status[train])

## Train and predict
cox = coxph(Surv(time, status) ~ ., data = data[train,])
p_cox_lp = predict(cox, newdata = data[test,])
p_cox_surv = survfit(cox, newdata = data[test,])
p_ranger = predict(
  ranger(Surv(time, status) ~ ., data = data[train, ]),
  data = data[test, ]
)
p_gbm = predict(
  blackboost(Surv(time, status) ~ ., data = data[train, ], family = Cindex()),
  newdata = data[test, ]
)

# Define Antolini's C-index
#  Copied from pycox
#   https://github.com/havakv/pycox/blob/master/pycox/evaluation/concordance.py
is_comparable = function(t_i, t_j, d_i, d_j) 
  ((t_i < t_j) & d_i) | ((t_i == t_j) & (d_i | d_j))

is_concordant = function(s_i, s_j, t_i, t_j, d_i, d_j)
  (s_i < s_j) & is_comparable(t_i, t_j, d_i, d_j)

sum_comparable = function(t, d) {
  count = 0
  for (i in seq_along(t)) {
    for (j in seq_along(t)) {
      if (j != i) {
        count = count + is_comparable(t[i], t[j], d[i], d[j])
      }
    }
  }
  count
}

sum_concordant_disc = function(s, t, d, s_idx) {
  count = 0
  for (i in seq_along(t)) {
    idx = s_idx[i]
    for (j in seq_along(t)) {
      if (j != i) {
        count = count +
          is_concordant(s[idx, i], s[idx, j], t[i], t[j], d[i], d[j])
      }
    }
  }
  count
}

# truth - Surv object corresponding to true survival outcomes
# surv - predicted survival matrix (T x n)
# surv_idx - 'surv_idx[i]' gives index in 'surv' corresponding to
#   the event time of individual 'i'.
antolini = function (truth, surv, surv_idx) {
  durations = truth[, "time"]
  events = truth[, "status"]
  sum_concordant_disc(surv, durations, events, surv_idx) /
    sum_comparable(durations, events)
}

## Calculative Cindex for CPH and GBM
# Harrell
harrell_cph = concordance(target ~ p_cox_lp, reverse = TRUE)$concordance
harrell_rsf = NA
harrell_gbm = concordance(target ~ p_gbm)$concordance

# Uno
uno_cph = concordance(target ~ p_cox_lp, reverse = TRUE, timewt = "n/G2")$concordance
uno_rsf = NA
uno_gbm = concordance(target ~ p_gbm, timewt = "n/G2")$concordance


## Method 1 - Ensemble mortality - higher value = more deaths = higher risk
ensemble_rsf = concordance(target ~ rowSums(p_ranger$chf), reverse = TRUE)$concordance
ensemble_cph = concordance(target ~ rowSums(-log(t(p_cox_surv$surv))), reverse = TRUE)$concordance
ensemble_gbm = NA

# Ensemble mortality for rsf
# This is based on the cumulative hazard function (CHF) of the RSF model
# which is calculated at 26 time-points
# These correspond to the unique non-censored times in the training set
dim(p_ranger$chf)
CHF_rsf <- p_ranger$chf
EM_rsf <- rowSums(CHF_rsf)

# Expected mortality for cox ph
# This is calculated based on the predicted survival function
# which is calculated at 51 time-points
# These correspond to the unique times (censored and non-censored) in the training set
dim(p_cox_surv$surv)
CHF_cox <- -log(t(p_cox_surv$surv))
EM_cox <- rowSums(CHF_cox)

# What if I restrict to only the CHF matching the 26 time-points as in RSF
cph.times <- sort(unique(data[train,]$time))
EM_cox_2 <- rowSums(CHF_cox[,which(cph.times %in% p_ranger$unique.death.times)])
EM_cox_3 <- rowSums(CHF_cox[,-which(cph.times %in% p_ranger$unique.death.times)])

plot(EM_cox_2 / EM_cox)
plot(EM_cox_3 / EM_cox)

plot((EM_cox_2 + EM_cox_3) / EM_cox)

# Nelson-Aalen estimator CHF

fit <- survfit(Surv(time, status) ~ 1, data = data[test,])
naest1 <- cumsum(fit$n.event/fit$n.risk)

## Method 2 - Antolini
antolini_rsf = antolini(
  target, t(p_ranger$survival),
  findInterval(target[, "time"], p_ranger$unique.death.times, all.inside = TRUE)
)
antolini_cph = antolini(
  target, p_cox_surv$surv,
  findInterval(target[, "time"], p_cox_surv$time, all.inside = TRUE)
)
antolini_gbm = NA


## Method 3 - Distribution summary (no extrapolation)
cox_surv = t(p_cox_surv$surv)
colnames(cox_surv) = p_cox_surv$time
ranger_surv = p_ranger$survival
colnames(ranger_surv) = p_ranger$unique.death.times

# Higher value = longer expected lifetime = lower risk - Absurd value due to improper distribution
summary_naive_cph = concordance(target ~ distr6::as.Distribution(1 - cox_surv, fun = "cdf")$mean())$concordance
summary_naive_rsf = concordance(target ~ distr6::as.Distribution(1 - ranger_surv, fun = "cdf")$mean())$concordance
summary_naive_gbm = NA

## Method 4 - Distribution summary (extrapolation)
cox_surv = cbind(1, cox_surv, 0) # Add probabilities 1 and 0
colnames(cox_surv)[1] = "0"
colnames(cox_surv)[ncol(cox_surv)] = tail(p_cox_surv$time, 1) + 1e-3
ranger_surv = cbind(1, ranger_surv, 0) # Add probabilities 1 and 0
colnames(ranger_surv)[1] = "0"
colnames(ranger_surv)[ncol(ranger_surv)] = tail(p_ranger$unique.death.times, 1) + 1e-3

summary_extr_cph = concordance(target ~ distr6::as.Distribution(1 - cox_surv, fun = "cdf")$mean())$concordance
summary_extr_rsf = concordance(target ~ distr6::as.Distribution(1 - ranger_surv, fun = "cdf")$mean())$concordance
summary_extr_gbm = NA

## Method 5 - Single probability comparison
distr_cox = distr6::as.Distribution(1 - cox_surv, fun = "cdf")
distr_rsf = distr6::as.Distribution(1 - ranger_surv, fun = "cdf")

# survival - higher value = higher prob survival = lower risk
cox_prob_concordance = rsf_prob_concordance = numeric(nrow(p_cox_surv$surv))
for (i in seq_along(cox_prob_concordance)) {
  cox_prob_concordance[i] = concordance(target ~ p_cox_surv$surv[i, ])$concordance
}
for (i in seq_along(rsf_prob_concordance)) {
  rsf_prob_concordance[i] = concordance(target ~ p_ranger$survival[, i])$concordance
}

rsf_max = which.max(rsf_prob_concordance)
rsf_min = which.min(rsf_prob_concordance)
rsf_rand = sample(seq_along(rsf_prob_concordance), 1)

prob_cph = cox_prob_concordance[c(rsf_min, rsf_max, rsf_rand)]
prob_rsf = rsf_prob_concordance[c(rsf_min, rsf_max, rsf_rand)]
prob_gbm = rep(NA, 3)

matrix(c(
  harrell_cph, harrell_rsf, harrell_gbm,
  uno_cph, uno_rsf, uno_gbm,
  antolini_cph, antolini_rsf, antolini_gbm,
  prob_cph[1], prob_rsf[1], prob_gbm[1],
  prob_cph[2], prob_rsf[2], prob_gbm[2],
  prob_cph[3], prob_rsf[3], prob_gbm[3],
  summary_naive_cph, summary_naive_rsf, summary_naive_gbm,
  summary_extr_cph, summary_extr_rsf, summary_extr_gbm,
  ensemble_cph, ensemble_rsf, ensemble_gbm
),
ncol = 3, byrow = TRUE,
dimnames = list(c(
  "Harrell", "Uno", "Antolini", "Prob (min)", "Prob (max)",
  "Prob (rand)", "Summary (naive)", "Summary (extr)", "Ensemble"
), c("CPH", "RSF", "GBM"))
)