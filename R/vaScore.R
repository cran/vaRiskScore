#' VA CVD Risk Score (2021)
#'
#' Calculates the cardiovascular (CVD) risk score for women military service
#' members and veterans.
#'
#' @param age Patient age (years: 1-110)
#' @param race Patient race (1 = White, 2 = African American, 3 = Hispanic)
#' @param TC Total cholesterol (mg/dL: 0-400)
#' @param HDL HDL cholesterol (mg/dL: 0-200)
#' @param SBP Systolic blood pressure (mmHg: 0-300)
#' @param SBPTrt Patient is on a blood pressure medication (1 = Yes, 0 = No)
#' @param Smoke Current smoker (1 = Yes, 0 = No)
#' @param DM Diabetes (1 = Yes, 0 = No)
#' @param Depression Major Depression (1 = Yes, 0 = No)
#' @param verbose logical: should input (patient profile) be printed.
#'
#' @return Estimated 10-year CVD Risk for VA women military service
#' members and veterans.
#'
#' @export
#'
#' @examples
#' library(vaRiskScore)
#' vaScore(age = 50,
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
#'Jeon‐Slaughter, H., Chen, X., Tsai, S., Ramanan, B., & Ebrahimi, R. (2021).
#'Developing an internally validated veterans affairs women cardiovascular
#'disease risk score using Veterans Affairs National Electronic Health Records.
#'Journal of the American Heart Association, 10(5), e019217.
#' @author
#' Xiaofei Chen; Haekyung Jeon‐Slaughter

vaScore <-
  function(age = 50,
           race = 1,
           SBPTrt = 1,
           SBP = 150,
           TC = 100,
           HDL = 50,
           DM = 1,
           Smoke = 1,
           Depression = 1,
           verbose = TRUE) {
    # race:   1=White, 2=African American; 3=Hispanic.
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

    if (age <= 0 |
        age >= 110 |
        SBP <= 0 |
        SBP >= 300 | TC < 0 | TC >= 400 | HDL <= 0 | HDL >= 200
        |
        !race %in% c(1:3) |
        !SBPTrt %in% c(0, 1) |
        !DM %in% c(0, 1) |
        !Smoke %in% c(0, 1) | !Depression %in% c(0, 1)) {
      return('Please double check input (see manual for allowable values)!')
    }

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

    getRisk = round(getRisk * 100, 2)
    risk_score = paste0("The predicted 10-year ASCVD risk is ", getRisk, "%")

    if (verbose) {
      cat("++++++++++++++++++++++++++++++++++++++++++++++++++",
          "\n")
      cat("Input:", "\n")
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
      cat(risk_score)
      return(getRisk/100)
    } else{
      return(getRisk/100)
    }

  }
