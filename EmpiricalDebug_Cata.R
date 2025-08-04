calculate.pairs.ties <- function(time, delta, mx, eval.times,
                                 tiedpredIn = 1, tiedoutcomeIn = 1, 
                                 tiedmatchIn = 1) {
  
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
  
  case1_pec <- 0
  case1.ties_pec <- 0
  case2_pec <- 0
  case2.ties_pec <- 0
  case3_pec <- 0
  case4_pec <- 0
  case5_pec <- 0
  case5.ties_pec <- 0
  case6_pec <- 0
  case6.ties_pec <- 0
  case7_pec <- 0
  concordant1_pec <- 0
  concordant1.ties_pec <- 0
  concordant2_pec <- 0
  concordant2.ties_pec <- 0
  concordant5_pec <- 0
  concordant5.ties_pec <- 0
  concordant6_pec <- 0
  concordant6.ties_pec <- 0
  
  
  ## calculate weights
  tmp <- data.frame(Y=time, status=delta, pred = mx)
  tmp <- tmp[order(tmp$Y),]
  
  weight_i <- pec::ipcw(formula=Surv(Y,status)~1,
                        data=tmp,
                        method="marginal",
                        times=unique(tmp$Y),
                        subjectTimes=tmp$Y, 
                        what = "IPCW.subjectTimes")$IPCW.subjectTimes
  
  weight_j <- pec::ipcw(formula=Surv(Y,status)~1,
                        data=tmp,
                        method="marginal",
                        times=unique(tmp$Y),
                        subjectTimes=tmp$Y,
                        subjectTimesLag=0,
                        what = "IPCW.times")$IPCW.times
  
  tindex = match(tmp$Y, unique(tmp$Y))
  
  time <- tmp$Y
  delta = tmp$status 
  mx = tmp$pred
  
  for(i in 1:(n - 1)) {
    for(j in (i + 1):n) {
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
          # pec
          wi <- weight_i[i]
          wj <- weight_j[tindex[i]]
          
          ww <- wi * wj
          if((time[i] <= eval.times) & (wi > 0 && wj > 0)){
            if (mx[i] == mx[j]) {
              if (tiedpredIn == 1) {
                case1.ties_pec <- case1.ties_pec + 1 / (ww)
                concordant1.ties_pec <- concordant1.ties_pec + 0.5 / (ww)
              }
            } else {
              case1_pec <- case1_pec + 1 / (ww)
              if (mx[i] > mx[j]) {
                concordant1_pec <- concordant1_pec + 1 / (ww)
              }
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
          wi <- weight_i[j]
          wj <- weight_j[tindex[j]]
          
          ww <- wi * wj
          if((time[j] <= eval.times) & (wi > 0 && wj > 0)){
            if (mx[i] == mx[j]) {
              if (tiedpredIn == 1) {
                case1.ties_pec <- case1.ties_pec + 1 / (ww)
                concordant1.ties_pec <- concordant1.ties_pec + 0.5 / (ww)
              }
            } else {
              case1_pec <- case1_pec + 1 / (ww)
              if (mx[i] < mx[j]) {
                concordant1_pec <- concordant1_pec + 1 / (ww)
              }
            }
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
          # pec
          wi <- weight_i[i]
          wj <- weight_j[tindex[i]]
          
          ww <- wi * wj
          if ((time[i] <= eval.times) & (wi > 0 && wj > 0)) {
            if (mx[i] == mx[j]) {
              if (tiedpredIn == 1) {
                case2.ties_pec <- case2.ties_pec + 1 / ww
                concordant2.ties_pec <- concordant2.ties_pec + 0.5 / ww
              }
            } else {
              case2_pec <- case2_pec + 1 / ww
              if (mx[i] > mx[j]) {
                concordant2_pec <- concordant2_pec + 1 / ww
              }
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
          # pec
          wi <- weight_i[j]
          wj <- weight_j[tindex[j]]
          
          ww <- wi * wj
          if ((time[j] <= eval.times) & (wi > 0 && wj > 0)) {
            if (mx[i] == mx[j]) {
              if (tiedpredIn == 1) {
                case2.ties_pec <- case2.ties_pec + 1 / ww
                concordant2.ties_pec <- concordant2.ties_pec + 0.5 / ww
              }
            } else {
              case2_pec <- case2_pec + 1 / ww
              if (mx[i] < mx[j]) {
                concordant2_pec <- concordant2_pec + 1 / ww
              }
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
        
        # Case 5 - only SurvMetrics considers this
        if((time[i] == time[j]) & (delta[i] == 1) & (delta[j] == 1)) {
          if(mx[i] == mx[j]) {
            case5.ties <- case5.ties + 1
            concordant5.ties <- concordant5.ties + 1 # adds 1 as per SurvMetrics
          } else {
            case5 <- case5 + 1    
            concordant5 <- concordant5 + 0.5 # adds 0.5 as per SurvMetrics
          }
          wi <- weight_i[i] ## is this or j ?
          wj <- weight_j[tindex[i]] ## is this or j ?
          ww <- wi * wj
          if ((time[i] <= eval.times) & (wi > 0 && wj > 0)) {
            if (mx[i] == mx[j]) {
              if (tiedmatchIn == 1) {
                case5_pec <- case5_pec + 1 / ww
                concordant5_pec <- concordant5_pec + 1 / ww
              } else if (tiedpredIn == 1 && tiedoutcomeIn == 1) {
                case5_pec <- case5_pec + 1 / ww
                concordant5.ties_pec <- concordant5.ties_pec + 0.5 / ww
              }
            } else if (tiedoutcomeIn == 1) {
              case5_pec <- case5_pec + 1 / ww
              if (mx[i] > mx[j]) {
                concordant5_pec <- concordant5_pec + 1 / ww
              }
            }
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
          wi <- weight_i[i]
          wj <- weight_j[tindex[i]]
          
          ww <- wi * wj
          if ((time[i] <= eval.times) & (wi > 0 && wj > 0)) {
            if (mx[i] > mx[j]) {
              case6_pec <- case6_pec + 1 / ww
              concordant6_pec <- concordant6_pec + 1 / ww
            }
            if (mx[i] == mx[j] && tiedpredIn == 1) {
              case6_pec <- case6_pec + 1 / ww
              concordant6.ties_pec <- concordant6.ties_pec + 0.5 / ww
            }
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
          wi <- weight_i[j]
          wj <- weight_j[tindex[j]]
          
          ww <- wi * wj
          if ((time[j] <= eval.times) & (wi > 0 && wj > 0)) {
            if (mx[i] == mx[j] && tiedpredIn == 1) {
              case6_pec <- case6_pec + 1 / ww
              concordant6.ties_pec <- concordant6.ties_pec + 0.5 / ww
            }
            if (mx[i] < mx[j]) {
              case6_pec <- case6_pec + 1 / ww
              concordant6_pec <- concordant6_pec + 1 / ww
            }
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
       # pec
       case1_pec = case1_pec,
       case1.ties_pec = case1.ties_pec,
       case2_pec = case2_pec,
       case2.ties_pec = case2.ties_pec,
       case5_pec = case5_pec,
       case5.ties_pec = case5.ties_pec,
       case6_pec = case6_pec,
       case6.ties_pec = case6.ties_pec,
       concordant1_pec = concordant1_pec,
       concordant1.ties_pec = concordant1.ties_pec,
       concordant2_pec = concordant2_pec,
       concordant2.ties_pec = concordant2.ties_pec,
       concordant5_pec = concordant5_pec,
       concordant5.ties_pec = concordant5.ties_pec,
       concordant6_pec = concordant6_pec,
       concordant6.ties_pec = concordant6.ties_pec,
       
       total_concordant = concordant1 + concordant2 + concordant3 +
         concordant4 + concordant5 + concordant6 + concordant7 +
         concordant1.ties + concordant2.ties + 
         concordant5.ties + concordant6.ties,
       total_cases = case1 + case2 + case3 + case4 + case5 + case6 + case7 +
         case1.ties + case2.ties + case5.ties + case6.ties)
}

calculate.pairs.ties.cata <- function(time, delta, mx, eval.times) {
  
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
  
  case1_pec <- 0
  case1.ties_pec <- 0
  case2_pec <- 0
  case2.ties_pec <- 0
  case3_pec <- 0
  case4_pec <- 0
  case5_pec <- 0
  case5.ties_pec <- 0
  case6_pec <- 0
  case6.ties_pec <- 0
  case7_pec <- 0
  concordant1_pec <- 0
  concordant1.ties_pec <- 0
  concordant2_pec <- 0
  concordant2.ties_pec <- 0
  concordant5_pec <- 0
  concordant5.ties_pec <- 0
  concordant6_pec <- 0
  concordant6.ties_pec <- 0
  
  
  ## calculate weights
  tmp <- data.frame(Y=time, status=delta, pred = mx)
  tmp <- tmp[order(tmp$Y),]
  
  weight_i <- pec::ipcw(formula=Surv(Y,status)~1,
                        data=tmp,
                        method="marginal",
                        times=unique(tmp$Y),
                        subjectTimes=tmp$Y, 
                        what = "IPCW.subjectTimes")$IPCW.subjectTimes
  
  weight_j <- pec::ipcw(formula=Surv(Y,status)~1,
                        data=tmp,
                        method="marginal",
                        times=unique(tmp$Y),
                        subjectTimes=tmp$Y,
                        subjectTimesLag=0,
                        what = "IPCW.times")$IPCW.times
  
  tindex = match(tmp$Y, unique(tmp$Y))
  
  time <- tmp$Y
  delta = tmp$status 
  mx = tmp$pred
  
  for(i in 1:(n - 1)) {
    for(j in (i + 1):n) {
      if((i != j)) {
        
        # pec weights
        wi <- weight_i[i]
        wj <- weight_j[tindex[i]]
        
        ww <- wi * wj
        
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
          
          if((time[i] <= eval.times) & (wi > 0 && wj > 0)){
            if (mx[i] == mx[j]) {
              case1.ties_pec <- case1.ties_pec + 1 / (ww)
              concordant1.ties_pec <- concordant1.ties_pec + 0.5 / (ww)
            } else {
              case1_pec <- case1_pec + 1 / (ww)
              if (mx[i] > mx[j]) {
                concordant1_pec <- concordant1_pec + 1 / (ww)
              }
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
          
          if((time[j] <= eval.times) & (wi > 0 && wj > 0)){
            if (mx[i] == mx[j]) {
              case1.ties_pec <- case1.ties_pec + 1 / (ww)
              concordant1.ties_pec <- concordant1.ties_pec + 0.5 / (ww)
            } else {
              case1_pec <- case1_pec + 1 / (ww)
              if (mx[i] < mx[j]) {
                concordant1_pec <- concordant1_pec + 1 / (ww)
              }
            }
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
          
          if ((time[i] <= eval.times) & (wi > 0 && wj > 0)) {
            if (mx[i] == mx[j]) {
              case2.ties_pec <- case2.ties_pec + 1 / ww
              concordant2.ties_pec <- concordant2.ties_pec + 0.5 / ww
            } else {
              case2_pec <- case2_pec + 1 / ww
              if (mx[i] > mx[j]) {
                concordant2_pec <- concordant2_pec + 1 / ww
              }
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
          
          if ((time[j] <= eval.times) & (wi > 0 && wj > 0)) {
            if (mx[i] == mx[j]) {
              case2.ties_pec <- case2.ties_pec + 1 / ww
              concordant2.ties_pec <- concordant2.ties_pec + 0.5 / ww
            } else {
              case2_pec <- case2_pec + 1 / ww
              if (mx[i] < mx[j]) {
                concordant2_pec <- concordant2_pec + 1 / ww
              }
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
        
        # Case 5 - only SurvMetrics considers this
        if((time[i] == time[j]) & (delta[i] == 1) & (delta[j] == 1)) {
          if(mx[i] == mx[j]) {
            case5.ties <- case5.ties + 1
            concordant5.ties <- concordant5.ties + 1 # adds 1 as per SurvMetrics
          } else {
            case5 <- case5 + 1    
            concordant5 <- concordant5 + 0.5 # adds 0.5 as per SurvMetrics
          }
          
          if ((time[i] <= eval.times) & (wi > 0 && wj > 0)) {
            
            if (mx[i] == mx[j]) {
              case5.ties_pec <- case5.ties_pec + 1 / ww
              concordant5.ties_pec <- concordant5.ties_pec + 0.5 / ww
            } else {
              case5_pec <- case5_pec + 1 / ww
              if (mx[i] > mx[j]) {
                concordant5_pec <- concordant5_pec + 1 / ww
              }
            }
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
          
          if ((time[i] <= eval.times) & (wi > 0 && wj > 0)) {
            
            if (mx[i] == mx[j]) {
              case6.ties_pec <- case6.ties_pec + 1 / ww
              concordant6.ties_pec <- concordant6.ties_pec + 0.5 / ww
            } else {
              case6_pec <- case6_pec + 1 / ww
              if (mx[i] > mx[j]) {
                concordant6_pec <- concordant6_pec + 1 / ww
              }
            }
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
          
          if ((time[j] <= eval.times) & (wi > 0 && wj > 0)) {
            
            if (mx[i] == mx[j]) {
              case6.ties_pec <- case6.ties_pec + 1 / ww
              concordant6.ties_pec <- concordant6.ties_pec + 0.5 / ww
            } else {
              case6_pec <- case6_pec + 1 / ww
              if (mx[i] < mx[j]) {
                concordant6_pec <- concordant6_pec + 1 / ww
              }
            }
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
       # pec
       case1_pec = case1_pec,
       case1.ties_pec = case1.ties_pec,
       case2_pec = case2_pec,
       case2.ties_pec = case2.ties_pec,
       case5_pec = case5_pec,
       case5.ties_pec = case5.ties_pec,
       case6_pec = case6_pec,
       case6.ties_pec = case6.ties_pec,
       concordant1_pec = concordant1_pec,
       concordant1.ties_pec = concordant1.ties_pec,
       concordant2_pec = concordant2_pec,
       concordant2.ties_pec = concordant2.ties_pec,
       concordant5_pec = concordant5_pec,
       concordant5.ties_pec = concordant5.ties_pec,
       concordant6_pec = concordant6_pec,
       concordant6.ties_pec = concordant6.ties_pec,
       
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

comparablelifelines <- function(MyC) {
  MyC$case1 + MyC$case1.ties +
    MyC$case2 + MyC$case2.ties +
    #MyC$case5 + MyC$case5.ties +
    MyC$case6 + MyC$case6.ties
}

concordantlifelines <- function(MyC) {
  MyC$concordant1 + MyC$concordant1.ties +
    MyC$concordant2 + MyC$concordant2.ties +
    #MyC$concordant5 + MyC$concordant5.ties +
    MyC$concordant6 + MyC$concordant6.ties 
}

comparablesksurv <- function(MyC) {
  MyC$case1 + MyC$case1.ties +
    MyC$case2 + MyC$case2.ties +
    #MyC$case5 + MyC$case5.ties +
    MyC$case6 + MyC$case6.ties
}

concordantsksurv <- function(MyC) {
  MyC$concordant1 + MyC$concordant1.ties +
    MyC$concordant2 + MyC$concordant2.ties +
    #MyC$concordant5 + MyC$concordant5.ties +
    MyC$concordant6 + MyC$concordant6.ties 
}


comparablepec <- function(MyC, 
                          tiedpredIn = 1, tiedoutcomeIn = 1, tiedmatchIn = 1) {
  
  out <- MyC$case1_pec + MyC$case2_pec
  
  if(tiedpredIn == 1) {
    out <- out + MyC$case1.ties_pec + MyC$case2.ties_pec
  } 
  
  if(tiedoutcomeIn == 0) {
    out <- out + MyC$case6_pec
    if(tiedmatchIn == 1 & tiedpredIn == 1) {
      out <- out + MyC$case5.ties_pec + MyC$case6.ties_pec
    } else if(tiedmatchIn == 1 & tiedpredIn == 0) {
      out <- out + MyC$case5.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 1) {
      out <- out + MyC$case6.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 0) {
      out <- out + 0
    }
  } else if(tiedoutcomeIn == 1) {
    out <- out + MyC$case5_pec + MyC$case6_pec
    if(tiedmatchIn == 1 & tiedpredIn == 1) {
      out <- out + MyC$case5.ties_pec + MyC$case6.ties_pec
    } else if(tiedmatchIn == 1 & tiedpredIn == 0) {
      out <- out + MyC$case5.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 1) {
      out <- out + MyC$case5.ties_pec + MyC$case6.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 0) {
      out <- out + 0
    }
  }  
  
  return(out)
}

concordantpec <- function(MyC, 
                          tiedpredIn = 1, tiedoutcomeIn = 1, tiedmatchIn = 1) {
  
  out <- MyC$concordant1_pec + MyC$concordant2_pec
  
  if(tiedpredIn == 1) {
    out <- out + MyC$concordant1.ties_pec + MyC$concordant2.ties_pec
  } 
  
  if(tiedoutcomeIn == 0) {
    out <- out + MyC$concordant6_pec
    if(tiedmatchIn == 1 & tiedpredIn == 1) {
      out <- out + 2*MyC$concordant5.ties_pec + MyC$concordant6.ties_pec
    } else if(tiedmatchIn == 1 & tiedpredIn == 0) {
      out <- out + 2*MyC$concordant5.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 1) {
      out <- out + MyC$concordant6.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 0) {
      out <- out + 0
    }
  } else if(tiedoutcomeIn == 1) {
    out <- out + MyC$concordant5_pec + MyC$concordant6_pec
    if(tiedmatchIn == 1 & tiedpredIn == 1) {
      out <- out + 2*MyC$concordant5.ties_pec + MyC$concordant6.ties_pec
    } else if(tiedmatchIn == 1 & tiedpredIn == 0) {
      out <- out + 2*MyC$concordant5.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 1) {
      out <- out + MyC$concordant5.ties_pec + MyC$concordant6.ties_pec
    } else if(tiedmatchIn == 0 & tiedpredIn == 0) {
      out <- out + 0
    }
  } 
  
  return(out)
}



# Setting 6: all possible types of ties

pec_loop <- function(time, death, age, MyC) {
  
  results <- expand.grid(
    tiedoutcomeIn = c(0, 1),
    tiedpredIn = c(0, 1),
    tiedmatchIn = c(0, 1)
  )
  results[, c("Comparable_pec", "Comparable_ours")] <- NA
  results[, c("Concordant_pec", "Concordant_ours")] <- NA
  results[, c("Cindex_pec", "Cindex_ours")] <- NA
  
  
  for(i in seq_len(nrow(results))) {
    # pec 
    pec_true <- calculatepec(time, death, age, 
                             tiedpredIn = results$tiedpredIn[i], 
                             tiedoutcomeIn = results$tiedoutcomeIn[i],
                             tiedmatchIn = results$tiedmatchIn[i])
    
    results$Comparable_pec[i] <- pec_true$Pairs$res / 2
    results$Concordant_pec[i] <- pec_true$Concordant$res / 2
    results$Cindex_pec[i] <- pec_true$AppCindex$res
    
    results$Comparable_ours[i] <- comparablepec(MyC, 
                                                tiedpredIn = results$tiedpredIn[i], 
                                                tiedoutcomeIn = results$tiedoutcomeIn[i],
                                                tiedmatchIn = results$tiedmatchIn[i])
    results$Concordant_ours[i] <- concordantpec(MyC,
                                                tiedpredIn = results$tiedpredIn[i], 
                                                tiedoutcomeIn = results$tiedoutcomeIn[i],
                                                tiedmatchIn = results$tiedmatchIn[i])
    results$Cindex_ours[i] <- results$Concordant_ours[i] / results$Comparable_ours[i]
    
  }
  
  results$Comparable_diff <- round(results$Comparable_pec - results$Comparable_ours, 1)
  results$Concordant_diff <- round(results$Concordant_pec - results$Concordant_ours, 1)
  results$Cindex_diff <- results$Cindex_pec - results$Cindex_ours
  
  return(results)
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

# Lifelines wrapper for R
lifelinesR <- function(time, predicted, censoring){
  # Load python libraries
  np <- import("numpy", convert = FALSE)
  lifelines_utils <- import("lifelines.utils", convert = FALSE)
  concordance_internal <- import("lifelines.utils.concordance", convert = FALSE)
  
  time_np      <- np$array(time, dtype = "float64")
  predicted_np <- np$array(1 - predicted, dtype = "float64")
  censoring_np <- np$array(as.integer(censoring), dtype = "bool")
  
  c_index <- lifelines_utils$concordance_index(
    event_times = time_np,
    predicted_scores = predicted_np,
    event_observed = censoring_np
  )
  
  stats <- concordance_internal$`_concordance_summary_statistics`(
    event_times = time_np,
    predicted_event_times = predicted_np,
    event_observed = censoring_np
  )

  list(
    Cindex     = py_to_r(c_index),
    concordant = py_to_r(stats[[0L]]),
    tied       = py_to_r(stats[[1L]]),
    comparable = py_to_r(stats[[2L]])
  )
}

sksurv.censoredR <- function(time, predicted, censoring, tied_tol = 1e-8) {
  np <- import("numpy", convert = FALSE)
  sksurv_internal <- import("sksurv.metrics", convert = FALSE)  # underscore

  time_py <- np$array(as.numeric(time), dtype = "float64")
  predicted_py <- np$array(as.numeric(predicted), dtype = "float64")
  censoring_py <- np$array(as.logical(censoring), dtype = "bool")
  weights_py <- np$ones_like(time_py, dtype = "float64") 
  
  result <- sksurv_internal$`_estimate_concordance_index`(
    event_indicator = censoring_py,
    event_time = time_py,
    estimate = predicted_py,
    weights = weights_py, ### the sksurv.ipcw should just change this
    tied_tol = tied_tol
  )
  
  list(
    Cindex      = py_to_r(result[[0L]]),
    concordant  = py_to_r(result[[1L]]),
    discordant  = py_to_r(result[[2L]]),
    tied_risk   = py_to_r(result[[3L]]),
    last_tied_time = py_to_r(result[[4L]])
  )
  
}

# Wrapper function to avoid code duplication
calculatepec <- function(time, death, age, 
                         tiedpredIn = 1, tiedoutcomeIn = 1, tiedmatchIn = 1) {
  tmp <- data.frame(age, death, time)
  eval.times <- max(time)
  pec::cindex(
    list("res" = as.matrix(1 - age)),
    formula = Hist(time, death) ~ 1, 
    data = tmp, 
    eval.times = eval.times,
    tiedPredictionsIn = I(tiedpredIn == 1),
    tiedOutcomeIn = I(tiedoutcomeIn == 1),
    tiedMatchIn = I(tiedmatchIn == 1),
    cens.model = "marginal", # default marginal
  )
}



cindexSRC_R <- function(Y, status, times, pred,
                        tiedpredIn = 1, tiedoutcomeIn = 1, tiedmatchIn = 1) {
  # Parameters:
  # Y: Vector of survival times
  # status: Event indicator (1 = event, 0 = censored)
  # times: Evaluation times
  # weight_i: IPCW weights for subjects
  # weight_j: IPCW weights for time points or subject pairs
  # pred: Matrix of predictions (subjects x times)
  # tindex: Time indices for survival times
  # tiedpredIn, tiedoutcomeIn, tiedmatchIn: Tie-handling parameters
  
  N <- length(Y)  # Number of subjects
  NT <- length(times)  # Number of evaluation times
  #status <- ifelse(status == 1, 0, 1)
  tmp <- data.frame(Y, status, pred)
  tmp <- tmp[order(tmp$Y),]

  weight_i <- pec::ipcw(formula=Surv(Y,status)~1,
                        data=tmp,
                        method="marginal",
                        times=unique(tmp$Y),
                        subjectTimes=tmp$Y, 
                        what = "IPCW.subjectTimes")$IPCW.subjectTimes
  
  #tmp$ipcw_subject_times <- weight.i
  
  weight_j <- pec::ipcw(formula=Surv(Y,status)~1,
                        data=tmp,
                        method="marginal",
                        times=unique(tmp$Y),
                        subjectTimes=tmp$Y,
                        subjectTimesLag=0,
                        what = "IPCW.times")$IPCW.times
  
  tindex = match(tmp$Y, unique(tmp$Y))
  #tmp$ipcw_times <- weight.j
  
  Y = tmp$Y
  status = tmp$status
  pred = as.matrix(1 - tmp$pred)
  
  # Initialize result vectors
  C <- numeric(NT)
  conc <- numeric(NT)
  pairs <- numeric(NT)
  
  # Loop over evaluation times
  for (s in seq_len(NT)) {
    conc[s] <- 0
    pairs[s] <- 0
    
    # Loop over subjects
    for (i in seq_len(N)) {
      # Usable pairs: i's event time must be <= current time and uncensored
      if (Y[i] <= times[s] && status[i] == 1) {
        for (j in seq(i + 1, N)) {  # Only consider pairs (i, j) with j > i
          # Compute IPCW weights
          wi <- weight_i[i]
          wj <- weight_j[tindex[i]]
          
          ww <- wi * wj
          
          # Skip pairs with zero weights
          if (wi > 0 && wj > 0) {
            # Tied outcomes and predictions: Fully concordant if tiedmatchIn == TRUE
            if (tiedmatchIn == 1 && Y[i] == Y[j] && status[j] == 1 &&
                pred[i, s] == pred[j, s]) {
              pairs[s] <- pairs[s] + 1 / ww
              conc[s] <- conc[s] + 1 / ww
            } else {
              # If tiedoutcomeIn == 0, exclude pairs with tied outcomes (unless j is censored)
              if (tiedoutcomeIn == 1 || (Y[i] != Y[j] || status[j] == 0)) {
                # Tied predictions: Include as half-concordant if tiedpredIn == TRUE
                if (pred[i, s] == pred[j, s]) {
                  if (tiedpredIn == 1) {
                    pairs[s] <- pairs[s] + 1 / ww
                    conc[s] <- conc[s] + 1 / (2 * ww)
                  }
                } else {
                  # Concordant if pred[i, s] < pred[j, s]
                  pairs[s] <- pairs[s] + 1 / ww
                  if (pred[i, s] < pred[j, s]) {
                    conc[s] <- conc[s] + 1 / ww
                  }
                }
              }
            }
          }
        }
      }
    }
    
    # Compute C-index for time s
    if (pairs[s] > 0) {
      C[s] <- conc[s] / pairs[s]
    } else {
      C[s] <- NA  # Avoid division by zero
    }
  }
  
  list(Cindex = C,
       concordant = conc, 
       comparable = pairs
  )
}




library(Hmisc)
library(SurvMetrics)
library(survival)
library(pec)

####### Setting 1: no ties ###### 
set.seed(1)
n <- 1000
age <- rnorm(n, 50, 10)
d.time <- rweibull(n, shape = 2, scale = 10)
cens   <- runif(n, 5, 20)
death  <- ifelse(d.time <= cens, 1, 0)
d.time <- pmin(d.time, cens)

MyC <- calculate.pairs.ties.cata(d.time, death, age, eval.times = max(d.time))

## Hmisc - OK
C1 <- Hmisc::rcorr.cens(-age, Surv(d.time, death))
comparableHmisc(MyC)
C1["Relevant Pairs"]/2# OK
concordantHmisc(MyC)
C1["Concordant"]/2 # OK
uncertainHmisc(MyC)
C1["Uncertain"]/2 # OK

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time, death), predicted = -age) # OK
C2$Concordant + C2$ConcordantTies # This is half of the concordant pairs in Hmisc
concordantSurvMetrics(MyC) 
C2$Comparable   # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC)
# This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time, death), predicted = -age)
concordantSurvMetrics(MyC) / comparableSurvMetrics(MyC) # OK

## Survival - OK
C3 <- survival::concordance(Surv(d.time, death) ~ age, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] # OK
concordantSurvival(MyC)  # OK
C3$count[["concordant"]] + C3$count[["discordant"]] 
comparableSurvival(MyC)
concordantSurvival(MyC)/comparableSurvival(MyC)

## Lifelines, we the same as in Survival 
C4 <- lifelinesR(d.time, age, death)
concordantlifelines(MyC)
C4$concordant + C4$tied*0.5
comparablelifelines(MyC)
C4$comparable
(C4$concordant + C4$tied*0.5) / C4$comparable
C4$Cindex

## sksurv
C5 <- sksurv.censoredR(d.time, age, death)
C5$Cindex
C5$concordant / (C5$concordant + C5$discordant)
C5$concordant
C5$concordant + C5$discordant
#C5$tied_risk
concordantsksurv(MyC)
comparablesksurv(MyC)
#C5$last_tied_time

## pec
pec_loop(d.time, death, age, MyC)

####### Setting 2: ties in time (only in uncensored) #####
d.time.ties1 <- round(d.time, 1)
death.ties1  <- d.time.ties1 <= cens
d.time.ties1 <- pmin(d.time.ties1, cens)

MyC <- calculate.pairs.ties.cata(d.time.ties1, death.ties1, age, eval.times = max(d.time.ties1))

# Hmisc - OK
C1 <- rcorr.cens(-age, Surv(d.time.ties1, death.ties1))
comparableHmisc(MyC)
C1["Relevant Pairs"]/2 # OK
concordantHmisc(MyC)
C1["Concordant"]/2 # OK
uncertainHmisc(MyC)
C1["Uncertain"]/2 # OK

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time.ties1, death.ties1), predicted = -age) # OK
C2$Concordant + C2$ConcordantTies # This is half of the concordant pairs in Hmisc
concordantSurvMetrics(MyC) 
C2$Comparable   # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) 
# This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time.ties1, death.ties1), predicted = -age)
concordantSurvMetrics(MyC) / comparableSurvMetrics(MyC) # OK

## Survival - OK but need to remove case 5 from the table
C3 <- survival::concordance(Surv(d.time.ties1, death.ties1) ~ age, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] # OK
concordantSurvival(MyC)  
C3$count[["concordant"]] + C3$count[["discordant"]] # OK
comparableSurvival(MyC)
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

## Lifelines, same results as survival
C4 <- lifelinesR(d.time.ties1, age, death.ties1)
concordantlifelines(MyC)
comparablelifelines(MyC)
C4$concordant + C4$tied*0.5
C4$comparable
(C4$concordant + C4$tied*0.5) / C4$comparable
C4$Cindex

## Sksurv
C5 <- sksurv.censoredR(d.time.ties1, age, death.ties1)
C5$Cindex
concordantsksurv(MyC) / comparablesksurv(MyC)
C5$concordant / (C5$concordant + C5$discordant)
concordantsksurv(MyC)
C5$concordant
comparablesksurv(MyC)
C5$concordant + C5$discordant
C5$tied_risk


## pec: 
pec_loop(d.time.ties1, death.ties1, age, MyC)


# Setting 3: ties in time (cens or uncensored)
d.time.ties2 <- round(d.time, 1)

MyC <- calculate.pairs.ties.cata(d.time.ties2, death, age,  eval.times = max(d.time.ties2))

## Hmisc - OK
C1 <- rcorr.cens(-age, Surv(d.time.ties2, death))
comparableHmisc(MyC)
C1["Relevant Pairs"] / 2 # OK
concordantHmisc(MyC)
C1["Concordant"] / 2# OK
uncertainHmisc(MyC)
C1["Uncertain"] / 2 # OK

## SurvMetrics - OK if we count case 6B as partially concordant (which I think doesn't make sense)
C2 <- Cindex_aux(Surv(d.time.ties2, death), predicted = -age)
C2$Concordant + C2$ConcordantTies  # OK
concordantSurvMetrics(MyC)+ MyC$concordant6.survmetrics  # OK
C2$Comparable   # OK
comparableSurvMetrics(MyC) # OK
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
concordantSurvival(MyC) 
C3$count[["concordant"]] + C3$count[["discordant"]] # OK
comparableSurvival(MyC) 
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

## Lifelines
C4 <- lifelinesR(d.time.ties2, age, death)
concordantlifelines(MyC)
C4$concordant
comparablelifelines(MyC)
C4$comparable
(C4$concordant + C4$tied*0.5) / C4$comparable
C4$Cindex

## sksurv
C5 <- sksurv.censoredR(d.time.ties2, age, death)
C5$Cindex
(concordantsksurv(MyC)) / (comparablesksurv(MyC))
C5$concordant / (C5$concordant + C5$discordant)
concordantsksurv(MyC)
C5$concordant
comparablesksurv(MyC)
C5$concordant + C5$discordant
C5$tied_risk


## pec: 
pec_loop(d.time.ties2, death, age, MyC)

# Setting 4: ties in predictions only
age_round <- round(age, 1)

MyC <- calculate.pairs.ties.cata(d.time, death, age_round, eval.times = max(d.time))

# Hmisc - OK
C1 <- rcorr.cens(-age_round, Surv(d.time, death))
comparableHmisc(MyC)
C1["Relevant Pairs"]/2 # OK
concordantHmisc(MyC)
C1["Concordant"]/2 # OK
uncertainHmisc(MyC)
C1["Uncertain"]/2 # OK

# Hmisc - exclude ties
C1 <- rcorr.cens(-age_round, Surv(d.time, death), outx = TRUE)
comparableHmisc(MyC) - MyC$case1.ties - MyC$case2.ties - MyC$case6.ties
C1["Relevant Pairs"]/2 # OK
concordantHmisc(MyC) - MyC$concordant1.ties - MyC$concordant2.ties - MyC$concordant6.ties
C1["Concordant"]/2 # OK
uncertainHmisc(MyC)
C1["Uncertain"]/2 # Not OK, but irrelevant for the Cindex (I just need to calculate the diff)

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time, death), predicted = -age_round) 
C2$Concordant + C2$ConcordantTies # This is half of the comparable pairs in Hmisc
concordantSurvMetrics(MyC) 
C2$Comparable # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) # This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time, death), predicted = -age_round)
concordantSurvMetrics(MyC) / comparableSurvMetrics(MyC) # OK

## Survival - OK
C3 <- survival::concordance(Surv(d.time, death) ~ age_round, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] + C3$count[["tied.x"]] / 2 # OK
concordantSurvival(MyC)
C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
comparableSurvival(MyC) 
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

# In survival, they first calculate somer's D and then transform into C-index
npair <- C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
somer <- (C3$count[["concordant"]] - C3$count[["discordant"]])/ npair
(somer +1)/2

## Lifelines
C4 <- lifelinesR(d.time, age_round, death)
concordantlifelines(MyC)
C4$concordant + C4$tied*0.5
comparablelifelines(MyC)
C4$comparable
(C4$concordant + C4$tied*0.5) / C4$comparable
C4$Cindex
C4$tied

# (num_correct + num_tied / 2) / num_pairs
# it is the same as survival
## sksurv
C5 <- sksurv.censoredR(d.time, age_round, death)
C5$Cindex
(concordantsksurv(MyC)/2) / (comparablesksurv(MyC)/2)
(C5$concordant + C5$tied_risk * 0.5) / (C5$concordant + C5$discordant + C5$tied_risk)
concordantsksurv(MyC)
(C5$concordant + C5$tied_risk * 0.5)
comparablesksurv(MyC)
C5$concordant + C5$discordant + C5$tied_risk
C5$tied_risk

## pec: 
pec_loop(d.time, death, age_round, MyC)


# Setting 5: ties in times or predictions (but not both)
MyC <- calculate.pairs.ties.cata(d.time.ties2, death, age_round, eval.times = max(d.time.ties2))

# Hmisc - OK
C1 <- rcorr.cens(-age_round, Surv(d.time.ties2, death))
comparableHmisc(MyC)
C1["Relevant Pairs"] /2# OK
concordantHmisc(MyC)
C1["Concordant"] /2# OK
uncertainHmisc(MyC)
C1["Uncertain"] /2# OK

## SurvMetrics - OK if I include case 6B as partially concordant (which I think doesn't make sense)
C2 <- Cindex_aux(Surv(d.time.ties2, death), predicted = -age_round)
C2$Concordant + C2$ConcordantTies 
# This only matches if I expand the partially concordant pairs (which doesn't make sense)
concordantSurvMetrics(MyC)  + MyC$concordant6.survmetrics 
C2$Comparable # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC)  # This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time.ties2, death), predicted = -age_round)
(concordantSurvMetrics(MyC) + MyC$concordant6.survmetrics) / comparableSurvMetrics(MyC) # OK

## Survival 
C3 <- survival::concordance(Surv(d.time.ties2, death) ~ age_round, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] + C3$count[["tied.x"]] / 2 # OK
concordantSurvival(MyC) 
C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
comparableSurvival(MyC) 
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

## Lifelines
C4 <- lifelinesR(d.time.ties2, age_round, death)
concordantlifelines(MyC)
(C4$concordant + C4$tied*0.5)
comparablelifelines(MyC)
C4$comparable
(C4$concordant + C4$tied*0.5) / C4$comparable
C4$Cindex

## sksurv
C5 <- sksurv.censoredR(d.time.ties2, age_round, death)
C5$Cindex
(concordantsksurv(MyC)) / (comparablesksurv(MyC))
(C5$concordant + C5$tied_risk * 0.5) / (C5$concordant + C5$discordant + C5$tied_risk)
concordantsksurv(MyC)
(C5$concordant + C5$tied_risk * 0.5)
comparablesksurv(MyC)
C5$concordant + C5$discordant + C5$tied_risk
C5$tied_risk



## pec: 
pec_loop(d.time.ties2, death, age_round, MyC)

# Setting 6: ties in times and predictions (simultaneously)
d.time.both <- c(d.time, rep(seq(0.1, 1.9, length.out = 10), each = 2))
death.both <- c(death, c(rep(c(1, 1), times = 3), rep(c(1, 0), times = 3), 
                         rep(c(0, 1), times = 2), rep(c(0, 0), times = 2)))
age.both <- c(age, rep(seq(30, 70, length.out = 10), each = 2))

MyC <- calculate.pairs.ties.cata(d.time.both, death.both, age.both, eval.times = max(d.time.both))

# Hmisc - OK
C1 <- rcorr.cens(-age.both, Surv(d.time.both, death.both))
comparableHmisc(MyC)
C1["Relevant Pairs"]/2 # OK
concordantHmisc(MyC)
C1["Concordant"] /2# OK
uncertainHmisc(MyC)
C1["Uncertain"] /2# OK

## SurvMetrics - OK
C2 <- Cindex_aux(Surv(d.time.both, death.both), predicted = -age.both)
C2$Concordant + C2$ConcordantTies # 
concordantSurvMetrics(MyC)+ MyC$concordant6.survmetrics 
C2$Comparable # This is half of the comparable pairs in Hmisc
comparableSurvMetrics(MyC) # This is half of the comparable ties in Hmisc
SurvMetrics::Cindex(Surv(d.time.both, death.both), predicted = -age.both)
(concordantSurvMetrics(MyC) + MyC$concordant6.survmetrics) / comparableSurvMetrics(MyC) # OK

## Survival 
C3 <- survival::concordance(Surv(d.time.both, death.both) ~ age.both, timewt = "n", reverse = TRUE)
C3$count[["concordant"]] + C3$count[["tied.x"]] /2# OK
concordantSurvival(MyC) 
C3$count[["concordant"]] + C3$count[["discordant"]] + C3$count[["tied.x"]]
comparableSurvival(MyC)
C3$concordance # OK
concordantSurvival(MyC) / comparableSurvival(MyC) 

## Lifelines
C4 <- lifelinesR(d.time.both,  age.both, death.both)
concordantlifelines(MyC)
(C4$concordant + C4$tied*0.5)
comparablelifelines(MyC)
C4$comparable
(C4$concordant + C4$tied*0.5) / C4$comparable
C4$Cindex

## sksurv
C5 <- sksurv.censoredR(d.time.both,  age.both, death.both)
C5$Cindex
(concordantsksurv(MyC)) / (comparablesksurv(MyC))
(C5$concordant + C5$tied_risk * 0.5) / (C5$concordant + C5$discordant + C5$tied_risk)
concordantsksurv(MyC)
(C5$concordant + C5$tied_risk * 0.5)
comparablesksurv(MyC)
C5$concordant + C5$discordant + C5$tied_risk
C5$tied_risk

## pec: 
pec_loop(d.time.both, death.both, age.both, MyC)


