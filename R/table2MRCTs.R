#' @name table2MRCTs
#'
#' @title Tables for two MRCTs paper
#'
#' @param numSimu The number of simulation times. Defaults to 100,000
#' @param choice Table number.
#' @param seed Random seed. Defaults to 1.
#'
#' @returns Tables shown in the reference
#'
#' @examples
#' table2MRCTs(100,1,1)
#'
#' @export
#'
table2MRCTs <- function(numSimu,choice,seed = 1) {
  #Table1
  table1 <- function(numSimu){
    consistencyProb <- data.frame()
    consistencyProbRow <- function(power,d,pCtrl){
      randRatio <- 1
      pi <- 0.5
      consistencyProbability <- 0.8
      alpha <- 0.025
      sigmaTrt <- sqrt((pCtrl+d)*(1-(pCtrl+d)))
      sigmaCtrl <- sqrt(pCtrl*(1-pCtrl))
      ctrlN <- (sigmaTrt^2/randRatio+sigmaCtrl^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
      ctrlN <- ceiling(ctrlN)
      trtN <- ceiling(randRatio*ctrlN)
      N <- trtN + ctrlN
      regionFraction <- regFrac(alpha = alpha, power = power, pi = pi, CP = consistencyProbability)
      count <- matrix(0,nrow = numSimu,ncol = 2)
      regionCtrlN <- ceiling(regionFraction*ctrlN)
      regionTrtN <- ceiling(regionFraction*trtN)
      actualRegionFraction <- regionCtrlN/ctrlN
      trtOutcome <- matrix(rbinom(trtN*numSimu,1,pCtrl+d),nrow = numSimu,ncol = trtN)
      ctrlOutcome <- matrix(rbinom(ctrlN*numSimu,1,pCtrl),nrow = numSimu,ncol = trtN)
      regionTrtOutcome <- trtOutcome[,1:regionTrtN]
      regionCtrlOutcome <- ctrlOutcome[,1:regionCtrlN]
      D <- apply(trtOutcome,1,mean) - apply(ctrlOutcome,1,mean)
      varD <- apply(trtOutcome,1,sd)^2/trtN + apply(ctrlOutcome,1,sd)^2/ctrlN
      z <- D/sqrt(varD)
      count[,1] <- ifelse(z > qnorm(1-alpha),1,0)
      D1 <- apply(regionTrtOutcome,1,mean) - apply(regionCtrlOutcome,1,mean)
      count[,2] <- ifelse(D1 > pi*D & z > qnorm(1-alpha),1,0)
      row <- c(power, d, pCtrl, N, format(regionFraction,digits = 3), format(sum(count[,2])/sum(count[,1]),digits = 3), ceiling(actualRegionFraction*1000)/1000)
      return(row)
    }
    for (power in c(0.8,0.9)) {
      for (d in c(0.1,0.15)) {
        for (pCtrl in c(0.5,0.6,0.7,0.8)) {
          consistencyProb <- rbind(consistencyProb, consistencyProbRow(power,d,pCtrl))
        }
      }
      for (d in c(0.2)) {
        for (pCtrl in c(0.5,0.6,0.7)) {
          consistencyProb <- rbind(consistencyProb, consistencyProbRow(power,d,pCtrl))
        }
      }
    }
    colnames(consistencyProb) <- c("power","d", "pCtrl", "N","f_k","empirical CP","actual f_k")
    return(consistencyProb)
  }
  #Table2
  table2 <- function(numSimu){
    consistencyProb <- data.frame()
    consistencyProbRow <- function(power,d){
      sigmaTrt <- 4
      sigmaCtrl <- sigmaTrt
      randRatio <- 1
      pi <- 0.5
      consistencyProbability <- 0.8
      alpha <- 0.025
      ctrlN <- (sigmaTrt^2/randRatio+sigmaCtrl^2)*(qnorm(1-alpha)+qnorm(power))^2/(d^2)
      ctrlN <- ceiling(ctrlN)
      trtN <- ceiling(randRatio*ctrlN)
      N <- trtN + ctrlN
      regionFraction <- regFrac(alpha = alpha, power = power, pi = pi, CP = consistencyProbability)
      count <- matrix(0,nrow = numSimu,ncol = 2)
      regionCtrlN <- ceiling(regionFraction*ctrlN)
      regionTrtN <- ceiling(regionFraction*trtN)
      actualRegionFraction <- regionCtrlN/ctrlN
      trtOutcome <- matrix(rnorm(trtN*numSimu,mean = d,sd = sigmaTrt),nrow = numSimu,ncol = trtN)
      ctrlOutcome <- matrix(rnorm(ctrlN*numSimu,mean = 0,sd = sigmaCtrl),nrow = numSimu,ncol = ctrlN)
      regionTrtOutcome <- trtOutcome[,1:regionTrtN]
      regionCtrlOutcome <- ctrlOutcome[,1:regionCtrlN]
      D <- apply(trtOutcome,1,mean) - apply(ctrlOutcome,1,mean)
      varD <- apply(trtOutcome,1,sd)^2/trtN + apply(ctrlOutcome,1,sd)^2/ctrlN
      z <- D/sqrt(varD)
      count[,1] <- ifelse(z > qnorm(1-alpha),1,0)
      D1 <- apply(regionTrtOutcome,1,mean) - apply(regionCtrlOutcome,1,mean)
      count[,2] <- ifelse(D1 > pi*D & z > qnorm(1-alpha),1,0)
      row <- c(power, d, N, format(regionFraction,digits = 3), format(sum(count[,2])/sum(count[,1]),digits = 3), ceiling(actualRegionFraction*1000)/1000)
      return(row)
    }
    for (power in c(0.8,0.9)) {
      for (d in c(1,1.5,2)) {
        consistencyProb <- rbind(consistencyProb, consistencyProbRow(power,d))
      }
    }
    colnames(consistencyProb) <- c("power","d", "N","f_k","empirical CP","actual f_k")
    return(consistencyProb)
  }
  #Table3
  table3 <- function(numSimu){
    consistencyProb <- data.frame()
    consistencyProbRow <- function(power1, power2, d, pCtrl1, pCtrl2){
      randRatio1 <- randRatio2 <- 1
      pi <- 0.5
      consistencyProbability <- 0.8
      alpha <- 0.025
      sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
      sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
      sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
      sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
      ctrlN1 <- (sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d^2)
      ctrlN1 <- ceiling(ctrlN1)
      trtN1 <- ceiling(randRatio1*ctrlN1)
      N1 <- trtN1 + ctrlN1
      ctrlN2 <- (sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d^2)
      ctrlN2 <- ceiling(ctrlN2)
      trtN2 <- ceiling(randRatio2*ctrlN2)
      N2 <- trtN2 + ctrlN2
      regionFraction <- regFrac2(alpha = alpha, power1 = power1, power2 = power2, pi = pi, CP = consistencyProbability)
      regionFraction1 <- regionFraction$rF1
      regionFraction2 <- regionFraction$rF2
      count <- matrix(0,nrow = numSimu,ncol = 2)
      regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
      regionTrtN1 <- ceiling(regionFraction1*trtN1)
      regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
      regionTrtN2 <- ceiling(regionFraction2*trtN2)
      w1 <- N1/(N1+N2)
      w2 <- 1 - w1
      actualRegionFraction1 <- regionCtrlN1/ctrlN1
      actualRegionFraction2 <- regionCtrlN2/ctrlN2
      trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
      ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
      trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
      ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
      regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
      regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
      regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
      regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
      D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
      D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
      varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
      varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
      z1 <- D1/sqrt(varD1)
      z2 <- D2/sqrt(varD2)
      count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
      D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
      Dpool <- w1*D1 + w2*D2
      Dpool1 <- w1*D11 + w2*D21
      count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      row <- c(power1,power2,d,pCtrl1, pCtrl2,N1,N2,format(regionFraction1,digits = 3),format(regionFraction2,digits = 3), format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
      # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
      return(row)
    }
    power1 <- 0.8
    for (power2 in c(0.8,0.9)) {
      for (d in c(0.1,0.15,0.2)) {
        pCtrl1 <- pCtrl2 <- 0.5
        consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, power2, d, pCtrl1, pCtrl2))
        pCtrl1 <- pCtrl2 <- 0.7
        consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, power2, d, pCtrl1, pCtrl2))
      }
    }
    power1 <- 0.9
    power2 <- 0.9
    for (d in c(0.1,0.15,0.2)) {
      pCtrl1 <- pCtrl2 <- 0.5
      consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, power2, d, pCtrl1, pCtrl2))
      pCtrl1 <- pCtrl2 <- 0.7
      consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, power2, d, pCtrl1, pCtrl2))
    }
    colnames(consistencyProb) <- c("power1","power2","d", "pCtrl1", "pCtrl2", "N1", "N2","f_k^1","f_k^2","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  #Table4
  table4 <- function(numSimu){
    consistencyProb <- data.frame()
    consistencyProbRow <- function(power1, power2, d){
      sigmaTrt1 <- 4
      sigmaCtrl1 <- sigmaTrt1
      sigmaTrt2 <- 4
      sigmaCtrl2 <- sigmaTrt2
      randRatio1 <- randRatio2 <- 1
      pi <- 0.5
      consistencyProbability <- 0.8
      alpha <- 0.025
      ctrlN1 <- (sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d^2)
      ctrlN1 <- ceiling(ctrlN1)
      trtN1 <- ceiling(randRatio1*ctrlN1)
      N1 <- trtN1 + ctrlN1
      ctrlN2 <- (sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d^2)
      ctrlN2 <- ceiling(ctrlN2)
      trtN2 <- ceiling(randRatio2*ctrlN2)
      N2 <- trtN2 + ctrlN2
      regionFraction <- regFrac2(alpha = alpha, power1 = power1, power2 = power2, pi = pi, CP = consistencyProbability)
      regionFraction1 <- regionFraction$rF1
      regionFraction2 <- regionFraction$rF2
      count <- matrix(0,nrow = numSimu,ncol = 2)
      regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
      regionTrtN1 <- ceiling(regionFraction1*trtN1)
      regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
      regionTrtN2 <- ceiling(regionFraction2*trtN2)
      w1 <- N1/(N1+N2)
      w2 <- 1 - w1
      actualRegionFraction1 <- regionCtrlN1/ctrlN1
      actualRegionFraction2 <- regionCtrlN2/ctrlN2
      trtOutcome1 <- matrix(rnorm(trtN1*numSimu,mean = d,sd = sigmaTrt1),nrow = numSimu,ncol = trtN1)
      ctrlOutcome1 <- matrix(rnorm(ctrlN1*numSimu,mean = 0,sd = sigmaCtrl1),nrow = numSimu,ncol = ctrlN1)
      trtOutcome2 <- matrix(rnorm(trtN2*numSimu,mean = d,sd = sigmaTrt2),nrow = numSimu,ncol = trtN2)
      ctrlOutcome2 <- matrix(rnorm(ctrlN2*numSimu,mean = 0,sd = sigmaCtrl2),nrow = numSimu,ncol = ctrlN2)
      regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
      regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
      regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
      regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
      D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
      D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
      varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
      varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
      z1 <- D1/sqrt(varD1)
      z2 <- D2/sqrt(varD2)
      count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
      D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
      Dpool <- w1*D1 + w2*D2
      Dpool1 <- w1*D11 + w2*D21
      count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      row <- c(power1,power2,d,N1,N2,format(regionFraction1,digits = 3),format(regionFraction2,digits = 3), format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
      # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
      return(row)
    }
    power1 <- 0.8
    for (power2 in c(0.8,0.9)) {
      for (d in c(1,1.5,2)) {
        consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, power2, d))
      }
    }
    power1 <- 0.9
    power2 <- 0.9
    for (d in c(1,1.5,2)) {
      consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, power2, d))
    }
    colnames(consistencyProb) <- c("power1","power2","d", "N1", "N2","f_k^1","f_k^2","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  #Table5
  table5 <- function(numSimu){
    consistencyProb <- data.frame()
    consistencyProbRow <- function(power1, randRatio1, power2, randRatio2, d, pCtrl1, pCtrl2){
      pi <- 0.5
      consistencyProbability <- 0.8
      alpha <- 0.025
      sigmaTrt1 <- sqrt((pCtrl1+d)*(1-(pCtrl1+d)))
      sigmaCtrl1 <- sqrt(pCtrl1*(1-pCtrl1))
      sigmaTrt2 <- sqrt((pCtrl2+d)*(1-(pCtrl2+d)))
      sigmaCtrl2 <- sqrt(pCtrl2*(1-pCtrl2))
      ctrlN1 <- (sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d^2)
      ctrlN1 <- ceiling(ctrlN1)
      trtN1 <- ceiling(randRatio1*ctrlN1)
      N1 <- trtN1 + ctrlN1
      ctrlN2 <- (sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d^2)
      ctrlN2 <- ceiling(ctrlN2)
      trtN2 <- ceiling(randRatio2*ctrlN2)
      N2 <- trtN2 + ctrlN2
      regionFraction <- regFrac2(alpha = alpha, power1 = power1, power2 = power2, pi = pi, d1 = d, d2 = d, sigmaTrt1 = sigmaTrt1, sigmaCtrl1 = sigmaCtrl1, sigmaTrt2 = sigmaTrt2, sigmaCtrl2 = sigmaCtrl2, randRatio1 = randRatio1, randRatio2 = randRatio2, CP = consistencyProbability)
      regionFraction1 <- regionFraction$rF1
      regionFraction2 <- regionFraction$rF2
      count <- matrix(0,nrow = numSimu,ncol = 2)
      regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
      regionTrtN1 <- ceiling(regionFraction1*trtN1)
      regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
      regionTrtN2 <- ceiling(regionFraction2*trtN2)
      w1 <- N1/(N1+N2)
      w2 <- 1 - w1
      actualRegionFraction1 <- regionCtrlN1/ctrlN1
      actualRegionFraction2 <- regionCtrlN2/ctrlN2
      trtOutcome1 <- matrix(rbinom(trtN1*numSimu,1,pCtrl1+d),nrow = numSimu,ncol = trtN1)
      ctrlOutcome1 <- matrix(rbinom(ctrlN1*numSimu,1,pCtrl1),nrow = numSimu,ncol = ctrlN1)
      trtOutcome2 <- matrix(rbinom(trtN2*numSimu,1,pCtrl2+d),nrow = numSimu,ncol = trtN2)
      ctrlOutcome2 <- matrix(rbinom(ctrlN2*numSimu,1,pCtrl2),nrow = numSimu,ncol = ctrlN2)
      regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
      regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
      regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
      regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
      D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
      D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
      varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
      varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
      z1 <- D1/sqrt(varD1)
      z2 <- D2/sqrt(varD2)
      count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
      D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
      Dpool <- w1*D1 + w2*D2
      Dpool1 <- w1*D11 + w2*D21
      count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      row <- c(power1,randRatio1,power2,randRatio2,d,pCtrl1, pCtrl2,N1,N2,format(regionFraction1,digits = 3),format(regionFraction2,digits = 3), format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
      # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
      return(row)
    }
    power1 <- 0.8
    randRatio1 <- 1
    power2 <- 0.8
    randRatio2 <- 2
    for (d in c(0.1,0.15,0.2)) {
      for (pCtrl1 in c(0.5,0.7)) {
        pCtrl2 <- pCtrl1
        consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, randRatio1, power2, randRatio2, d, pCtrl1, pCtrl2))
      }
    }
    power1 <- 0.8
    randRatio1 <- 1
    power2 <- 0.9
    randRatio2 <- 2
    for (d in c(0.1,0.15,0.2)) {
      for (pCtrl1 in c(0.5,0.7)) {
        pCtrl2 <- pCtrl1
        consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, randRatio1, power2, randRatio2, d, pCtrl1, pCtrl2))
      }
    }
    power1 <- 0.9
    randRatio1 <- 1
    power2 <- 0.9
    randRatio2 <- 2
    for (d in c(0.1,0.15,0.2)) {
      for (pCtrl1 in c(0.5,0.7)) {
        pCtrl2 <- pCtrl1
        consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, randRatio1, power2, randRatio2, d, pCtrl1, pCtrl2))
      }
    }
    colnames(consistencyProb) <- c("power1","randRatio1","power2","randRatio2","d", "pCtrl1", "pCtrl2", "N1", "N2","f_k^1","f_k^2","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  #Table6
  table6 <- function(numSimu){
    consistencyProb <- data.frame()
    consistencyProbRow <- function(power1, randRatio1, power2, randRatio2, d){
      sigmaTrt1 <- 4
      sigmaCtrl1 <- sigmaTrt1
      sigmaTrt2 <- 4
      sigmaCtrl2 <- sigmaTrt2
      pi <- 0.5
      consistencyProbability <- 0.8
      alpha <- 0.025
      ctrlN1 <- (sigmaTrt1^2/randRatio1+sigmaCtrl1^2)*(qnorm(1-alpha)+qnorm(power1))^2/(d^2)
      ctrlN1 <- ceiling(ctrlN1)
      trtN1 <- ceiling(randRatio1*ctrlN1)
      N1 <- trtN1 + ctrlN1
      ctrlN2 <- (sigmaTrt2^2/randRatio2+sigmaCtrl2^2)*(qnorm(1-alpha)+qnorm(power2))^2/(d^2)
      ctrlN2 <- ceiling(ctrlN2)
      trtN2 <- ceiling(randRatio2*ctrlN2)
      N2 <- trtN2 + ctrlN2
      regionFraction <- regFrac2(alpha = alpha, power1 = power1, power2 = power2, pi = pi, d1 = d, d2 = d, sigmaTrt1 = sigmaTrt1, sigmaCtrl1 = sigmaCtrl1, sigmaTrt2 = sigmaTrt2, sigmaCtrl2 = sigmaCtrl2, randRatio1 = randRatio1, randRatio2 = randRatio2, CP = consistencyProbability)
      regionFraction1 <- regionFraction$rF1
      regionFraction2 <- regionFraction$rF2
      count <- matrix(0,nrow = numSimu,ncol = 2)
      regionCtrlN1 <- ceiling(regionFraction1*ctrlN1)
      regionTrtN1 <- ceiling(regionFraction1*trtN1)
      regionCtrlN2 <- ceiling(regionFraction2*ctrlN2)
      regionTrtN2 <- ceiling(regionFraction2*trtN2)
      w1 <- N1/(N1+N2)
      w2 <- 1 - w1
      actualRegionFraction1 <- regionCtrlN1/ctrlN1
      actualRegionFraction2 <- regionCtrlN2/ctrlN2
      trtOutcome1 <- matrix(rnorm(trtN1*numSimu,mean = d,sd = sigmaTrt1),nrow = numSimu,ncol = trtN1)
      ctrlOutcome1 <- matrix(rnorm(ctrlN1*numSimu,mean = 0,sd = sigmaCtrl1),nrow = numSimu,ncol = ctrlN1)
      trtOutcome2 <- matrix(rnorm(trtN2*numSimu,mean = d,sd = sigmaTrt2),nrow = numSimu,ncol = trtN2)
      ctrlOutcome2 <- matrix(rnorm(ctrlN2*numSimu,mean = 0,sd = sigmaCtrl2),nrow = numSimu,ncol = ctrlN2)
      regionTrtOutcome1 <- trtOutcome1[,1:regionTrtN1]
      regionCtrlOutcome1 <- ctrlOutcome1[,1:regionCtrlN1]
      regionTrtOutcome2 <- trtOutcome2[,1:regionTrtN2]
      regionCtrlOutcome2 <- ctrlOutcome2[,1:regionCtrlN2]
      D1 <- apply(trtOutcome1,1,mean) - apply(ctrlOutcome1,1,mean)
      D2 <- apply(trtOutcome2,1,mean) - apply(ctrlOutcome2,1,mean)
      varD1 <- apply(trtOutcome1,1,sd)^2/trtN1 + apply(ctrlOutcome1,1,sd)^2/ctrlN1
      varD2 <- apply(trtOutcome2,1,sd)^2/trtN2 + apply(ctrlOutcome2,1,sd)^2/ctrlN2
      z1 <- D1/sqrt(varD1)
      z2 <- D2/sqrt(varD2)
      count[,1] <- ifelse(z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      D11 <-  apply(regionTrtOutcome1,1,mean) - apply(regionCtrlOutcome1,1,mean)
      D21 <-  apply(regionTrtOutcome2,1,mean) - apply(regionCtrlOutcome2,1,mean)
      Dpool <- w1*D1 + w2*D2
      Dpool1 <- w1*D11 + w2*D21
      count[,2] <- ifelse(Dpool1 > pi*Dpool & z1 > qnorm(1-alpha) & z2 > qnorm(1-alpha),1,0)
      row <- c(power1,randRatio1,power2,randRatio2,d,N1,N2,format(regionFraction1,digits = 3),format(regionFraction2,digits = 3), format(sum(count[,2])/sum(count[,1]),digits = 3),ceiling(actualRegionFraction1*1000)/1000,ceiling(actualRegionFraction2*1000)/1000)
      # consistencyProb equals the sum of the 1st col divided by the 2nd col, the last entry above
      return(row)
    }
    power1 <- 0.8
    randRatio1 <- 1
    power2 <- 0.8
    randRatio2 <- 2
    for (d in c(1,1.5,2)) {
      consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, randRatio1, power2, randRatio2, d))
    }
    power1 <- 0.8
    randRatio1 <- 1
    power2 <- 0.9
    randRatio2 <- 2
    for (d in c(1,1.5,2)) {
      consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, randRatio1, power2, randRatio2, d))
    }
    power1 <- 0.9
    randRatio1 <- 1
    power2 <- 0.9
    randRatio2 <- 2
    for (d in c(1,1.5,2)) {
      consistencyProb <- rbind(consistencyProb, consistencyProbRow(power1, randRatio1, power2, randRatio2, d))
    }
    colnames(consistencyProb) <- c("power1","randRatio1","power2","randRatio2","d", "N1", "N2","f_k^1","f_k^2","empirical CP","actual f_k^1","actual f_k^2")
    return(consistencyProb)
  }
  ####
  functionList <- list(table1,table2,table3,table4,table5,table6)
  if (choice >= 1 && choice <= length(functionList)) {
    set.seed(seed)
    functionList[[choice]](numSimu)
  } else {
    stop("Invalid choice!")
  }
}
