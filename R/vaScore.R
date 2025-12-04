#' VA CVD Risk Score (2021)
#'
#' Calculates the cardiovascular (CVD) risk score for civilian women, women military service
#' members and veterans. 
#'
#' @param population Military service (1 = Veteran, 2 = Current/Active, 3 = No Military)
#' @param age Patient age (years: 1-110)
#' @param race Patient race (1 = White, 2 = African American, 3 = Hispanic, 4 = Asian/Hawaiian Americans, 5 = American Indian/Indigenous Americans and Alaska Native)
#' @param TC Total cholesterol (mg/dL: 0-400)
#' @param HDL HDL cholesterol (mg/dL: 0-200)
#' @param SBP Systolic blood pressure (mmHg: 0-300)
#' @param SBPTrt Patient is on a blood pressure medication (1 = Yes, 0 = No)
#' @param Smoke Current smoker (1 = Yes, 0 = No)
#' @param DM Diabetes (1 = Yes, 0 = No)
#' @param Depression Major Depression (1 = Yes, 0 = No)
#' @param verbose logical: should input (patient profile) be printed.
#'
#' @return Estimated 10-year CVD Risk 
#'
#' @export
#'
#' @examples
#' library(vaRiskScore)
#' vaScore(population = 1,
#'         age = 50,
#'         race = 1,
#'         SBPTrt = 1,
#'         SBP = 150,
#'         TC = 100,
#'         HDL = 50,
#'         DM = 1,
#'         Smoke = 1,
#'         Depression = 1,
#'         verbose = TRUE)
#' @references
#'Jeon-Slaughter, H., Chen, X., Whyne, E.Z., Tsai, S., Barbosa, M.R., 
#'Ramanan, B., Bhushan, S., & Cao, D. (2025). External Validation of the 
#'Veterans Affairs Women Cardiovascular Disease Risk Score to
#'Nonveteran Women. JACC Advances, 4(9).
#' 
#'Jeon‐Slaughter, H., Chen, X., Tsai, S., Ramanan, B., & Ebrahimi, R. (2021).
#'Developing an internally validated veterans affairs women cardiovascular
#'disease risk score using Veterans Affairs National Electronic Health Records.
#'Journal of the American Heart Association, 10(5), e019217.
#' @author
#' Xiaofei Chen; Haekyung Jeon‐Slaughter

vaScore <-
  function(population = 1,
           age = 50,
           race = 1,
           SBPTrt = 1,
           SBP = 150,
           TC = 100,
           HDL = 50,
           DM = 1,
           Smoke = 1,
           Depression = 1,
           verbose = TRUE) {
    # population: Military service (1 = Veteran, 2 = Current/Active, 3 = No Military)
    # race: 1 = White, 2 = African American, 3 = Hispanic, 4 = Asian/Hawaiian Americans, 5 = American Indian/Indigenous Americans and Alaska Native
    # SBPTrt: 1=Yes;   2=No; same for others.

    x_age    =  age * 12
    x_lnage  =  log(x_age)
    x_lnagesq =  x_lnage ^ 2
    x_race   =  race
    x_SBPTrt =  SBPTrt
    x_SBP    =  SBP
    x_lnSBP  =  log(x_SBP)
    x_TC     =  TC
    x_lnTC   =  log(x_TC)
    x_HDL    =  HDL
    x_lnHDL  =  log(x_HDL)
    x_DM     =  DM
    x_SMK    =  Smoke
    x_Dep    =  Depression

    if (!population %in% c(1:3) | 
        age <= 0 |
        age >= 110 |
        SBP <= 0 |
        SBP >= 300 | TC < 0 | TC >= 400 | HDL <= 0 | HDL >= 200
        |
        !race %in% c(1:5) |
        !SBPTrt %in% c(0, 1) |
        !DM %in% c(0, 1) |
        !Smoke %in% c(0, 1) | !Depression %in% c(0, 1)) {
      return('Please double check input (see manual for allowable values)!')
    }
    if(population==3 & race %in% c(4:5)){
      cat("Note: current version does not support No Military Asian or American Indian","\n") 
      return('Please double check input (see manual for allowable values)!')
    }
#----------2021 publications for veteran White, AA, and Hispanic----------#
if(population==1 & race %in% c(1:3)){   
    if (x_race == 1) {
      Mlnage     = log(47.27 * 12)
      Mlnagesq   = Mlnage ^ 2
      MlnTotChol = log(198.63)
      MlnHDL     = log(53.48)
      MlnBP      = log(124.69)

      tmp <-
        exp(
          2.399 * x_lnage + 0.024 * x_lnTC - 1.350 * x_lnHDL - 0.208 * x_lnSBP * x_SBPTrt  +
            1.008 * x_lnSBP + 0.072 * x_SMK + 0.425 * x_DM + 0.244 *
            x_Dep +
            1.263 * x_SBPTrt - (
              2.399 * Mlnage + 0.024 * MlnTotChol - 1.350 * MlnHDL + 1.008 * MlnBP 
            )
        )
      getRisk    = 1 - 0.9438 ^ tmp
      raceCat    = "White"
    }


    if (x_race == 2) {
      Mlnage     = log(45.49 * 12)
      Mlnagesq   = Mlnage ^ 2
      MlnTotChol = log(192.09)
      MlnHDL     = log(56.69)
      MlnBP      = log(128.02)

      tmp <-
        exp(
          2.058 * x_lnage + 0.180 * x_lnTC - 1.339 * x_lnHDL + 1.246 * x_lnSBP * (x_SBPTrt) +
            0.411 * x_lnSBP - 0.020 * (x_SMK) + 0.276 * (x_DM) + 0.231 * (x_Dep) -
            5.795 * (x_SBPTrt) - (
              2.058 * Mlnage + 0.180 * MlnTotChol - 1.339 * MlnHDL + 0.411 * MlnBP 
            )
        )
      getRisk = 1 - 0.9442 ^ tmp
      raceCat = "African American"
    }


    if (x_race == 3) {
      Mlnage     = log(44.64 * 12)
      Mlnagesq   = Mlnage ^ 2
      MlnTotChol = log(195.47)
      MlnHDL     = log(53.85)
      MlnBP      = log(123.39)

      tmp <-
        exp(
          2.191 * x_lnage + 0.099 * x_lnTC - 1.225 * x_lnHDL - 3.714 * x_lnSBP * (x_SBPTrt)  +
            0.653 * x_lnSBP + 0.356 * (x_SMK) + 0.315 * (x_DM) + 0.311 * (x_Dep) +
            18.290 * (x_SBPTrt) - (
              2.191 * Mlnage + 0.099 * MlnTotChol - 1.225 * MlnHDL + 0.653 * MlnBP 
            )
        )
      getRisk = 1 - 0.9542 ^ tmp
      raceCat = "Hispanic"
    }
}
    
#------2025 JACC-ADV for active and civilian: White, AA, and Hispanic-------#    
if(population %in% c(2:3) & race %in% c(1:3)){   
    S10_pop2_W = 0.996
    S10_pop2_AA = 0.995
    S10_pop2_H = 0.997
    S10_pop3_W = 0.920
    S10_pop3_AA = 0.890
    S10_pop3_H = 0.907
    
      if (x_race == 1) {
        if(population==2){
        Mlnage     = log(32.37 * 12)
        Mlnagesq   = Mlnage ^ 2
        MlnTotChol = log(184.65)
        MlnHDL     = log(60.65)
        MlnBP      = log(118.75)
        }
        if(population==3){
          Mlnage     = log(45.90 * 12)
          Mlnagesq   = Mlnage ^ 2
          MlnTotChol = log(197.05)
          MlnHDL     = log(58.53)
          MlnBP      = log(119.24)
        }
        tmp <-
          exp(
            2.399 * x_lnage + 0.024 * x_lnTC - 1.350 * x_lnHDL - 0.208 * x_lnSBP * x_SBPTrt  +
              1.008 * x_lnSBP + 0.072 * x_SMK + 0.425 * x_DM + 0.244 *
              x_Dep +
              1.263 * x_SBPTrt - (
                2.399 * Mlnage + 0.024 * MlnTotChol - 1.350 * MlnHDL + 1.008 * MlnBP 
              )
          )
       if(population==2){getRisk    = 1 - S10_pop2_W ^ tmp}
       if(population==3){getRisk    = 1 - S10_pop3_W ^ tmp}
       raceCat = "White"
      }
      
      
      if (x_race == 2) {
        if(population==2){
          Mlnage     = log(33.46 * 12)
          Mlnagesq   = Mlnage ^ 2
          MlnTotChol = log(184.99)
          MlnHDL     = log(60.62)
          MlnBP      = log(121.15)
        }
        if(population==3){
          Mlnage     = log(43.60 * 12)
          Mlnagesq   = Mlnage ^ 2
          MlnTotChol = log(187.08)
          MlnHDL     = log(56.17)
          MlnBP      = log(127.92)
        }
        
        tmp <-
          exp(
            2.058 * x_lnage + 0.180 * x_lnTC - 1.339 * x_lnHDL + 1.246 * x_lnSBP * (x_SBPTrt) +
              0.411 * x_lnSBP - 0.020 * (x_SMK) + 0.276 * (x_DM) + 0.231 * (x_Dep) -
              5.795 * (x_SBPTrt) - (
                2.058 * Mlnage + 0.180 * MlnTotChol - 1.339 * MlnHDL + 0.411 * MlnBP 
              )
          )
        if(population==2){getRisk = 1 - S10_pop2_AA ^ tmp}
        if(population==3){getRisk = 1 - S10_pop3_AA ^ tmp}
        raceCat = "African American"
      }
      
      
      if (x_race == 3) {
        if(population==2){
          Mlnage     = log(27.48 * 12)
          Mlnagesq   = Mlnage ^ 2
          MlnTotChol = log(182.49)
          MlnHDL     = log(60.49)
          MlnBP      = log(115.84)
        }
        if(population==3){
          Mlnage     = log(39.63 * 12)
          Mlnagesq   = Mlnage ^ 2
          MlnTotChol = log(189.38)
          MlnHDL     = log(50.96)
          MlnBP      = log(116.84)
        }
        
        tmp <-
          exp(
            2.191 * x_lnage + 0.099 * x_lnTC - 1.225 * x_lnHDL - 3.714 * x_lnSBP * (x_SBPTrt)  +
              0.653 * x_lnSBP + 0.356 * (x_SMK) + 0.315 * (x_DM) + 0.311 * (x_Dep) +
              18.290 * (x_SBPTrt) - (
                2.191 * Mlnage + 0.099 * MlnTotChol - 1.225 * MlnHDL + 0.653 * MlnBP 
              )
          )
        if(population==2){getRisk = 1 - S10_pop2_H ^ tmp}
        if(population==3){getRisk = 1 - S10_pop3_H ^ tmp}
        raceCat = "Hispanic"
      }
}   
    
#------2025 JACC-ASIA for veteran and active: Asian and AI-------#    
if(population %in% c(1:2) & race %in% c(4:5)){   
      S10_pop1_Asian = 0.9729124
      S10_pop1_AI    = 0.9628521
      S10_pop2_Asian = 0.9729124
      S10_pop2_AI    = 0.9628521
      
      if (x_race == 4) {
        Mlnage     = log(34.91 * 12)
        Mlnagesq   = Mlnage ^ 2
        MlnTotChol = log(190.36)
        MlnHDL     = log(55.66)
        MlnBP      = log(119.35)
        
        tmp <-
          exp(
            2.662 * x_lnage + 0.313 * x_lnTC - 0.939 * x_lnHDL + 0.266 * x_lnSBP * x_SBPTrt  +
              0.368 * x_lnSBP -0.324 * x_SMK + 0.057 * x_DM + 0.220 *
              x_Dep - 0.434 * x_SBPTrt - (
                2.662 * Mlnage + 0.313 * MlnTotChol - 0.939 * MlnHDL + 0.368 * MlnBP 
              )
          )
        if(population==1){getRisk    = 1 - S10_pop1_Asian ^ tmp}
        if(population==2){getRisk    = 1 - S10_pop2_Asian ^ tmp}
        raceCat = "Asian/Hawaiian Americans"
      }
      
      
      if (x_race == 5) {
        Mlnage     = log(36.40 * 12)
        Mlnagesq   = Mlnage ^ 2
        MlnTotChol = log(190.90)
        MlnHDL     = log(53.23)
        MlnBP      = log(121.09)
        
        tmp <-
          exp(
            2.540 * x_lnage + 0.009 * x_lnTC - 0.560 * x_lnHDL - 2.247 * x_lnSBP * (x_SBPTrt) +
              1.311 * x_lnSBP + 0.127 * (x_SMK) + 0.266 * (x_DM) + 0.057 * (x_Dep) +
              1.174 * (x_SBPTrt) - (
                2.540 * Mlnage + 0.009 * MlnTotChol - 0.560 * MlnHDL + 1.311 * MlnBP 
              )
          )
        if(population==1){getRisk    = 1 - S10_pop1_AI ^ tmp}
        if(population==2){getRisk    = 1 - S10_pop2_AI ^ tmp}
        raceCat = "American Indian/Indigenous Americans and Alaska Native"
      }
      
}   
    
    if(population==3 & race %in% c(4:5)) getRisk = NA
    getRisk = round(getRisk * 100, 2)
    risk_score = paste0("The predicted 10-year ASCVD risk is ", getRisk, "%")

    
    if (verbose) {
      cat("++++++++++++++++++++++++++++++++++++++++++++++++++",
          "\n")
      cat("Input:", "\n")
      cat(paste("population: (1 = Veteran, 2 = Current/Active, 3 = No Military):", population), "\n")
      cat(paste("Age (year):", age), "\n")
      cat(paste("Race:", raceCat), "\n")
      cat(paste("Systolic blood pressure (mmHg):", SBP), "\n")
      cat(paste("Total cholesterol (mg/dL):", TC), "\n")
      cat(paste("HDL cholesterol (mg/dL):", HDL), "\n")
      cat(paste("Blood pressure treatment:", ifelse(SBPTrt == 1, 'Yes', 'No')), "\n")
      cat(paste("Diabetes:", ifelse(DM == 1, 'Yes', 'No')), "\n")
      cat(paste("Current smoker:", ifelse(Smoke == 1, 'Yes', 'No')), "\n")
      cat(paste("Major depression:", ifelse(Depression == 1, 'Yes', 'No')), "\n")
      cat("++++++++++++++++++++++++++++++++++++++++++++++++++",
          "\n")
      cat("Result:", "\n")
      cat(risk_score, "\n")
      return(getRisk/100)
    } else{
      return(getRisk/100)
    }

  }
