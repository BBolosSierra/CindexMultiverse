calculate.pairs.ties <- function(time, delta, mx) {
  
  n <- length(time)
  
  ## Total number of possible pairs
  total <- n*(n - 1)
  ## Calculate number of pairs under each case
  case1 <- 0
  case1.ties <- 0
  case2 <- 0
  case2.ties <- 0
  case3 <- 0
  case4 <- 0
  case5 <- 0
  case5.ties <- 0
  case6 <- 0
  case6.ties <- 0
  case7 <- 0
  concordant1 <- 0
  concordant1.ties <- 0
  concordant2 <- 0
  concordant2.ties <- 0
  concordant3 <- 0
  concordant4 <- 0
  concordant5 <- 0
  concordant5.ties <- 0
  concordant6 <- 0
  concordant6.ties <- 0
  concordant7 <- 0
  concordant6.survmetrics <- 0 # for SurvMetrics
  
  for(i in 1:n) {
    for(j in 1:n) {
      if((i != j)) {
        
        # Case 1 
        if((time[i] < time[j]) & (delta[i] == 1) & (delta[j] == 1))  {
          # We will separate comparable pairs such that predictions are
          # also tied as they are treated differently in diff methods
          if(mx[i] == mx[j]) {
            case1.ties <- case1.ties + 1
            concordant1.ties <- concordant1.ties + 0.5
          } else {
            case1 <- case1 + 1   
            if(mx[i] > mx[j]) {
              concordant1 <- concordant1 + 1   
            }
          }
        }
        if((time[i] > time[j]) & (delta[i] == 1) & (delta[j] == 1)) {
          if(mx[i] == mx[j]) {
            case1.ties <- case1.ties + 1
            concordant1.ties <- concordant1.ties + 0.5
          } else {
            case1 <- case1 + 1   
            if(mx[i] < mx[j])
              concordant1 <- concordant1 + 1  
          }
        }
        
        # Case 2
        if((time[i] < time[j]) & (delta[i] == 1) & (delta[j] == 0)) {
          if(mx[i] == mx[j]) {
            case2.ties <- case2.ties + 1
            concordant2.ties <- concordant2.ties + 0.5
          } else {
            case2 <- case2 + 1 
            if(mx[i] > mx[j]) {
              concordant2 <- concordant2 + 1  
            }
          }         
        }
        if((time[i] > time[j]) & (delta[i] == 0) & (delta[j] == 1)) {
          if(mx[i] == mx[j]) {
            case2.ties <- case2.ties + 1
            concordant2.ties <- concordant2.ties + 0.5
          } else {
            case2 <- case2 + 1 
            if(mx[i] < mx[j]) {
              concordant2 <- concordant2 + 1  
            }
          }
        }
        
        # Case 3
        if((time[i] < time[j]) & (delta[i] == 0) & (delta[j] == 1))  
          case3 <- case3 + 1
        if((time[i] > time[j]) & (delta[i] == 1) & (delta[j] == 0))  
          case3 <- case3 + 1
        
        # Case 4
        if((time[i] < time[j]) & (delta[i] == 0) & (delta[j] == 0))  
          case4 <- case4 + 1
        if((time[i] > time[j]) & (delta[i] == 0) & (delta[j] == 0))  
          case4 <- case4 + 1

        # Case 5 - only SurMetrics considers this
        if((time[i] == time[j]) & (delta[i] == 1) & (delta[j] == 1)) {
          if(mx[i] == mx[j]) {
            case5.ties <- case5.ties + 1
            concordant5.ties <- concordant5.ties + 1 # adds 1 as per SurvMetrics
          } else {
            case5 <- case5 + 1    
            concordant5 <- concordant5 + 0.5 # adds 0.5 as per SurvMetrics
          }
        }

        # Case 6
        if((time[i] == time[j]) & (delta[i] == 1) & (delta[j] == 0)) {
          if(mx[i] == mx[j]) {
            case6.ties <- case6.ties + 1
            concordant6.ties <- concordant6.ties + 0.5
          } else {
            case6 = case6 + 1   
            if(mx[i] > mx[j])
              concordant6 <- concordant6 + 1
            if(mx[i] < mx[j])
              concordant6.survmetrics <- concordant6.survmetrics + 0.5
          }
        }
        if((time[i] == time[j]) & (delta[i] == 0) & (delta[j] == 1)) {
          if(mx[i] == mx[j]) {
            case6.ties <- case6.ties + 1
            concordant6.ties <- concordant6.ties + 0.5
          } else {
            case6 <- case6 + 1   
            if(mx[i] < mx[j])
              concordant6 <- concordant6 + 1
            if(mx[i] > mx[j])
              concordant6.survmetrics <- concordant6.survmetrics + 0.5
          }
        }
        
        # Case 7
        if((time[i] == time[j]) & (delta[i] == 0) & (delta[j] == 0))  
          case7 <- case7 + 1       
      }  
    }
  }
  

  list(case1 = case1,
       case1.ties = case1.ties,
       case2 = case2,
       case2.ties = case2.ties,
       case3 = case3,
       case4 = case4,
       case5 = case5,
       case5.ties = case5.ties,
       case6 = case6,
       case6.ties = case6.ties,
       case7 = case7,
       concordant1 = concordant1,
       concordant1.ties = concordant1.ties,
       concordant2 = concordant2,
       concordant2.ties = concordant2.ties,
       concordant3 = concordant3,
       concordant4 = concordant4,
       concordant5 = concordant5,
       concordant5.ties = concordant5.ties,
       concordant6 = concordant6,
       concordant6.ties = concordant6.ties,
       concordant7 = concordant7,
       concordant6.survmetrics = concordant6.survmetrics,
       total_concordant = concordant1 + concordant2 + concordant3 +
         concordant4 + concordant5 + concordant6 + concordant7 +
         concordant1.ties + concordant2.ties + 
         concordant5.ties + concordant6.ties,
       total_cases = case1 + case2 + case3 + case4 + case5 + case6 + case7 +
         case1.ties + case2.ties + case5.ties + case6.ties)
}

comparableHmisc <- function(MyC) {
  MyC$case1 + MyC$case1.ties +
    MyC$case2 + MyC$case2.ties +
    MyC$case6 + MyC$case6.ties
}

concordantHmisc <- function(MyC) {
  MyC$concordant1 + MyC$concordant1.ties +
    MyC$concordant2 + MyC$concordant2.ties +
    MyC$concordant6 + MyC$concordant6.ties
}

uncertainHmisc <- function(MyC) {
  MyC$case3 + MyC$case4 + MyC$case7
}

comparableSurvMetrics <- function(MyC) {
  MyC$case1 + MyC$case1.ties +
    MyC$case2 + MyC$case2.ties +
    MyC$case5 + MyC$case5.ties +
    MyC$case6 + MyC$case6.ties
}

concordantSurvMetrics <- function(MyC) {
  MyC$concordant1 + MyC$concordant1.ties +
    MyC$concordant2 + MyC$concordant2.ties +
    MyC$concordant5 + MyC$concordant5.ties +
    MyC$concordant6 + MyC$concordant6.ties
}

comparableSurvival <- function(MyC) {
  MyC$case1 + MyC$case1.ties +
    MyC$case2 + MyC$case2.ties +
    #MyC$case5 + MyC$case5.ties +
    MyC$case6 + MyC$case6.ties
}

concordantSurvival <- function(MyC) {
  MyC$concordant1 + MyC$concordant1.ties +
    MyC$concordant2 + MyC$concordant2.ties +
    #MyC$concordant5 + MyC$concordant5.ties +
    MyC$concordant6 + MyC$concordant6.ties 
}

## Adaptation of SurvMetrics function to break the calculation into cases
Cindex_aux <- function(object, predicted) 
{
  time <- object[, 1]
  status <- object[, 2]
  if (length(time) != length(status)) {
    stop("The lengths of time and status are not equal")
  }
  if (length(time) != length(predicted)) {
    stop("The lengths of time and predicted are not equal")
  }
  if (missing(time) | missing(status) | missing(predicted)) {
    stop("Input value is missing")
  }
  if (any(is.na(time) | is.na(status) | is.na(predicted))) {
    stop("The input vector cannot have NA")
  }
  permissible <- 0
  concord <- 0
  concord12 <- 0
  concord5 <- 0
  concord6 <- 0
  par_concord <- 0
  par_concord12 <- 0
  par_concord5 <- 0
  par_concord6 <- 0
  n <- length(time)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      # Case 3 and 4
      if ((time[i] < time[j] & status[i] == 0) | (time[j] < 
                                                  time[i] & status[j] == 0)) {
        next
      }
      # Case 7
      if (time[i] == time[j] & status[i] == 0 & status[j] == 
          0) {
        next
      }
      permissible <- permissible + 1
      if (time[i] != time[j]) {
        # Case 1 A/B + 2 A/B
        if ((time[i] < time[j] & predicted[i] < predicted[j]) | 
            (time[j] < time[i] & predicted[j] < predicted[i])) {
          concord <- concord + 1
          concord12 <- concord12 + 1
        }
        # Case 1C + 2C
        else if (predicted[i] == predicted[j]) {
          par_concord <- par_concord + 0.5
          par_concord12 <- par_concord12 + 0.5
        }
      }
      # Case 5
      if (time[i] == time[j] & status[i] == 1 & status[j] == 
          1) {
        if (predicted[i] == predicted[j]) {
          concord <- concord + 1
          concord5 <- concord5 + 1
        }
        else {
          par_concord <- par_concord + 0.5
          par_concord5 <- par_concord5 + 0.5
        }
      }
      # Case 6
      if (time[i] == time[j] & ((status[i] == 1 & status[j] == 
                                 0) | (status[i] == 0 & status[j] == 1))) {
        if ((status[i] == 1 & predicted[i] < predicted[j]) | 
            (status[j] == 1 & predicted[j] < predicted[i])) {
          concord <- concord + 1
          concord6 <- concord6 + 1
        }
        else {
          #cat(paste0("Pair ", i, " and ", j, ".\n"))
          par_concord <- par_concord + 0.5
          par_concord6 <- par_concord6 + 0.5
        }
      }
    }
  }
  C_index <- (concord + par_concord)/permissible
  names(C_index) <- "C index"
  list("Cindex" = C_index,
       "Concordant" = concord,
       "Comparable" = permissible,
       "ConcordantTies" = par_concord,
       "Concordant12" = concord12,
       "ConcordantTies12" = par_concord12,
       "Concordant5" = concord5,
       "ConcordantTies5" = par_concord5,
       "Concordant6" = concord6,
       "ConcordantTies6" = par_concord6
  )
}

library(Hmisc)
library(SurvMetrics)
library(survival)

# Setting 1: no ties
set.seed(1)
n <- 1000
age <- rnorm(n, 50, 10)
d.time <- rweibull(n, shape = 2, scale = 10)
cens   <- runif(n, 5, 20)
death  <- ifelse(d.time <= cens, 1, 0)
d.time <- pmin(d.time, cens)

MyC <- calculate.pairs.ties(d.time, death, age)

## Hmisc - OK
C1 <- Hmisc::rcorr.cens(-age, Surv(d.time, death))
comparableHmisc(MyC)
C1["Relevant Pairs"] # OK
concordantHmisc(MyC)
C1["Concordant"] # OK
uncertainHmisc(MyC)
C1["Uncertain"] # OK

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time, death), predicted = -age) # OK
C2$Concordant + C2$ConcordantTies # This is half of the concordant pairs in Hmisc
concordantSurvMetrics(MyC) / 2
C2$Comparable   # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) / 2
# This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time, death), predicted = -age)
concordantSurvMetrics(MyC) / comparableSurvMetrics(MyC) # OK

## Survival - OK
C3 <- survival::concordance(Surv(d.time, death) ~ age, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] # OK
concordantSurvival(MyC) / 2 # OK
C3$count[["concordant"]] + C3$count[["discordant"]] 
comparableSurvival(MyC) / 2
concordantSurvival(MyC)/comparableSurvival(MyC)

# Setting 2: ties in time (only in uncensored)
d.time.ties1 <- round(d.time, 1)
death.ties1  <- d.time.ties1 <= cens
d.time.ties1 <- pmin(d.time.ties1, cens)

MyC <- calculate.pairs.ties(d.time.ties1, death.ties1, age)

# Hmisc - OK
C1 <- rcorr.cens(-age, Surv(d.time.ties1, death.ties1))
comparableHmisc(MyC)
C1["Relevant Pairs"] # OK
concordantHmisc(MyC)
C1["Concordant"] # OK
uncertainHmisc(MyC)
C1["Uncertain"] # OK

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time.ties1, death.ties1), predicted = -age) # OK
C2$Concordant + C2$ConcordantTies # This is half of the concordant pairs in Hmisc
concordantSurvMetrics(MyC) / 2
C2$Comparable   # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) / 2
# This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time.ties1, death.ties1), predicted = -age)
concordantSurvMetrics(MyC) / comparableSurvMetrics(MyC) # OK

## Survival - OK but need to remove case 5 from the table
C3 <- survival::concordance(Surv(d.time.ties1, death.ties1) ~ age, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] # OK
concordantSurvival(MyC) / 2 
C3$count[["concordant"]] + C3$count[["discordant"]] # OK
comparableSurvival(MyC) / 2
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

# Setting 3: ties in time (cens or uncensored)
d.time.ties2 <- round(d.time, 1)

MyC <- calculate.pairs.ties(d.time.ties2, death, age)

## Hmisc - OK
C1 <- rcorr.cens(-age, Surv(d.time.ties2, death))
comparableHmisc(MyC)
C1["Relevant Pairs"] # OK
concordantHmisc(MyC)
C1["Concordant"] # OK
uncertainHmisc(MyC)
C1["Uncertain"] # OK

## SurvMetrics - OK if we count case 6B as partially concordant (which I think doesn't make sense)
C2 <- Cindex_aux(Surv(d.time.ties2, death), predicted = -age)
C2$Concordant + C2$ConcordantTies  # OK
concordantSurvMetrics(MyC) / 2 + MyC$concordant6.survmetrics / 2 # OK
C2$Comparable   # OK
comparableSurvMetrics(MyC) / 2 # OK
SurvMetrics::Cindex(Surv(d.time.ties2, death), predicted = -age)
(concordantSurvMetrics(MyC) + MyC$concordant6.survmetrics) / comparableSurvMetrics(MyC) # OK
# Check example pairs - this suggests there is a bug or a weird treatment of ties
mypair <- c(951, 991)
cbind(d.time.ties2, death, age)[mypair,]
mypair <- c(927, 947)
cbind(d.time.ties2, death, age)[mypair,]

## Survival - OK but need to remove case 5 from the table
C3 <- survival::concordance(Surv(d.time.ties2, death) ~ age, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] # OK
concordantSurvival(MyC) / 2 
C3$count[["concordant"]] + C3$count[["discordant"]] # OK
comparableSurvival(MyC) / 2
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

# In survival, they first calculate somer's D and then transform into C-index
npair <- C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.y"]]
somer <- (C3$count[["concordant"]] - C3$count[["discordant"]])/ npair
(somer +1)/2

# Setting 4: ties in predictions only
age_round <- round(age, 1)

MyC <- calculate.pairs.ties(d.time, death, age_round)

# Hmisc - OK
C1 <- rcorr.cens(-age_round, Surv(d.time, death))
comparableHmisc(MyC)
C1["Relevant Pairs"] # OK
concordantHmisc(MyC)
C1["Concordant"] # OK
uncertainHmisc(MyC)
C1["Uncertain"] # OK

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time, death), predicted = -age_round) 
C2$Concordant + C2$ConcordantTies # This is half of the comparable pairs in Hmisc
concordantSurvMetrics(MyC) / 2
C2$Comparable # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) / 2 # This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time, death), predicted = -age_round)
concordantSurvMetrics(MyC) / comparableSurvMetrics(MyC) # OK

## Survival - OK
C3 <- survival::concordance(Surv(d.time, death) ~ age_round, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] + C3$count[["tied.x"]] / 2 # OK
concordantSurvival(MyC) / 2 
C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
comparableSurvival(MyC) / 2
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

# In survival, they first calculate somer's D and then transform into C-index
npair <- C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
somer <- (C3$count[["concordant"]] - C3$count[["discordant"]])/ npair
(somer +1)/2

# Setting 5: ties in times or predictions (but not both)
MyC <- calculate.pairs.ties(d.time.ties2, death, age_round)

# Hmisc - OK
C1 <- rcorr.cens(-age_round, Surv(d.time.ties2, death))
comparableHmisc(MyC)
C1["Relevant Pairs"] # OK
concordantHmisc(MyC)
C1["Concordant"] # OK
uncertainHmisc(MyC)
C1["Uncertain"] # OK

## SurvMetrics - OK if I include case 6B as partially concordant (which I think doesn't make sense)
C2 <- Cindex_aux(Surv(d.time.ties2, death), predicted = -age_round)
C2$Concordant + C2$ConcordantTies 
# This only matches if I expand the partially concordant pairs (which doesn't make sense)
concordantSurvMetrics(MyC) / 2 + MyC$concordant6.survmetrics / 2
C2$Comparable # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) / 2 # This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time.ties2, death), predicted = -age_round)
(concordantSurvMetrics(MyC) + MyC$concordant6.survmetrics) / comparableSurvMetrics(MyC) # OK

## Survival 
C3 <- survival::concordance(Surv(d.time.ties2, death) ~ age_round, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] + C3$count[["tied.x"]] / 2 # OK
concordantSurvival(MyC) / 2 
C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
comparableSurvival(MyC) / 2
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

# Setting 6: ties in times and predictions (simultaneously)
d.time.both <- c(d.time, rep(seq(0.1, 1.9, length.out = 10), each = 2))
death.both <- c(death, c(rep(c(1, 1), times = 3), rep(c(1, 0), times = 3), 
                         rep(c(0, 1), times = 2), rep(c(0, 0), times = 2)))
age.both <- c(age, rep(seq(30, 70, length.out = 10), each = 2))

MyC <- calculate.pairs.ties(d.time.both, death.both, age.both)

# Hmisc - OK
C1 <- rcorr.cens(-age.both, Surv(d.time.both, death.both))
comparableHmisc(MyC)
C1["Relevant Pairs"] # OK
concordantHmisc(MyC)
C1["Concordant"] # OK
uncertainHmisc(MyC)
C1["Uncertain"] # OK

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time.both, death.both), predicted = -age.both)
C2$Concordant + C2$ConcordantTies # 
concordantSurvMetrics(MyC) / 2 + MyC$concordant6.survmetrics / 2
C2$Comparable # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) / 2 # This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time.both, death.both), predicted = -age.both)
(concordantSurvMetrics(MyC) + MyC$concordant6.survmetrics) / comparableSurvMetrics(MyC) # OK

## Survival 
C3 <- survival::concordance(Surv(d.time.both, death.both) ~ age.both, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] + C3$count[["tied.x"]] / 2 # OK
concordantSurvival(MyC) / 2 
C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
comparableSurvival(MyC) / 2
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 
