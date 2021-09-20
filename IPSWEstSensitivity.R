



#############################################################################################
# R code for paper: Estimating the Optimal Timing of Surgery From Observational Data        #
# Date:             2020-01-18                                                              #
# R version:        3.6.0 (2019-04-26) -- "Planting of a Tree"                              #
# OS:               macOS X Catalina                                                        #
# Platform:         x86_64-apple-darwin15.6.0 (64-bit)                                      #  
#############################################################################################

source('DesignMatrixFun.R')

# Load R packages

require(dplyr)
require(discSurv)
require(VGAM)
require(likert)
require(coxphw)
require(smoothHR)
require(sqldf)


#----------------------------------- public data set ---------------------------------------- # 

# Load data from NIH/NHLBI Pediatric Heart Network,  
# who provided the SVRT data that we used in the preparation of this work.
# https://www.pediatricheartnetwork.com/pud_login.asp?study_id=SVR
# We call the raw data set downloaded from above site 'raw'.

#----------------------------- create event time and indicator ----------------------------- # 

data <- raw

# Delete one patient with missing S2P info

data <- subset(data, blind_id != 94)

# Initialize global variables for Norwood and S2P (Z and T)

event1_time <- rep(NA, nrow(data)) 
event2_time <- rep(NA, nrow(data))
censor1     <- rep(NA, nrow(data))
censor2     <- rep(NA, nrow(data))


for (i in 1:nrow(data)) {
  if (data$HAVESTG2[i] == '1:Yes') {
    event1_time[i]   <- data$STG2_age[i] - data$norw_age[i] + 1
    censor1[i]       <- 0   # Censoring indicator: 0:Censoring; 1:Not censoring
    
    if (data$transplant[i] == '2:No' & data$death[i] == '2:No') {
      event2_time[i] <- data$svrend_age[i] - data$STG2_age[i] + 1
      censor2[i]     <- 0
    }
    
    if (data$transplant[i] == '2:No' & data$death[i] == '1:Yes') {
      event2_time[i] <- data$death_age[i] - data$STG2_age[i] + 1
      censor2[i]     <- 1
    }
    if (data$transplant[i] == '1:Yes') {
      event2_time[i] <- data$transplant_age[i] - data$STG2_age[i] + 1
      censor2[i]     <- 1
    }
  }
  
  if (data$HAVESTG2[i] == '2:No') {
    if (data$transplant[i] == '1:Yes') {
      event1_time[i] <- data$transplant_age[i] - data$norw_age[i] + 1
      censor1[i]     <- 1
    }
    if (data$transplant[i] == '2:No' & data$death[i] == '1:Yes') {
      event1_time[i] <- data$death_age[i] - data$norw_age[i] + 1
      censor1[i]     <- 1
    }
    event2_time[i]   <- NA
    censor2[i]       <- NA
  }
  
}

data <- cbind(data,event1_time,censor1,event2_time,censor2)

# Create Race2 category: 1:White; 2:American African; 0:others.

RACE2                                           <- rep(0,nrow(data))
RACE2[data$RACE=='1:White']                     <- 1
RACE2[data$RACE=='2:Black or African American'] <- 2
data                                            <- cbind(data, RACE2)
rm(RACE2)

# Select patients with S2P

data2                                           <- data[data$HAVESTG2=='1:Yes',]

#----------------------------- calculate generalized propensity score ----------------------------- # 

# Number of competing risks

N_CR         <- 3
eventCR_time <- as.numeric(cut(data$event1_time, breaks = c(seq(0,366,30.5), 700),labels =1:13))

# Initialize global variables for competing-risk indicators

censorCR     <- rep(0, nrow(data))
censorCR0    <- rep(0, nrow(data))
censorCR1    <- rep(0, nrow(data)) 
censorCR2    <- rep(0, nrow(data)) 
censorCR3    <- rep(0, nrow(data))

for(i in 1:nrow(data)) {
  if (data$HAVESTG2[i] == '1:Yes' & data$ELECTIVE[i] == '1:Yes') {
    censorCR[i] <- 2  # S2P
  } else if (data$transplant[i] == '1:Yes' &
             data$HAVESTG2[i] == '2:No' &  data$censor1[i] == 1)
  {
    censorCR[i] <- 3
  }  # transplant
  else if (data$transplant[i] == '2:No' &
           data$HAVESTG2[i] == '2:No' &  data$censor1[i] == 1)
  {
    censorCR[i] <- 1
  }  # death
}

for(i in 1:length(censorCR)){
  if(censorCR[i]==0){
    censorCR0[i] <- 1    # censorCR0: non-elective S2P
  }
  if(censorCR[i]==1){
    censorCR1[i] <- 1    # censorCR1: death
  }
  if(censorCR[i]==2){
    censorCR2[i] <- 1    # censorCR2: elective S2P
  }
  if(censorCR[i]==3){
    censorCR3[i] <- 1    # censorCR3: heart transplant
  }
}

# Sensitivity analysis with competing risks:
# elecitve-S2P, non-elective S2P, death/heart transplant 

censorCR4   <- censorCR1 + censorCR3

data        <- cbind(data, eventCR_time, censorCR0, censorCR2, censorCR4)

dataLong    <-
  dataLongCompRisks(
    dataSet      = data,
    timeColumn   = "eventCR_time",
    eventColumns = c('censorCR0','censorCR2','censorCR4')
  )

dataLong$eT <-  dataLong$e1 + dataLong$e2 + dataLong$e3

# Create S2P timing category: 

Z           <- rep(NA,nrow(dataLong))

for(i in 1:nrow(dataLong)) {
  if (as.numeric(dataLong[i, 'timeInt']) <  4)
    Z[i] <- 1
  if (as.numeric(dataLong[i, 'timeInt']) == 4)
    Z[i] <- 2
  if (as.numeric(dataLong[i, 'timeInt']) == 5)
    Z[i] <- 3
  if (as.numeric(dataLong[i, 'timeInt']) == 6)
    Z[i] <- 4
  if (as.numeric(dataLong[i, 'timeInt']) == 7)
    Z[i] <- 5
  if (as.numeric(dataLong[i, 'timeInt']) >  7)
    Z[i] <- 6
}
dataLong$Z   <- Z


# Fit the multinomial logit model for discrete-time competing risks by Tutz

multLogitVGM <-
  vgam(
    cbind(e0, e1, e2, e3) ~ timeInt + factor(GENDER) + factor(RACE2) + factor(ctrt) + BWT + reverse.levels(factor(PRENATDX)) +
      reverse.levels(factor(ATRESIA)) + reverse.levels(factor(OBSTRUCT)) + nordis_age,
    family = multinomial(refLevel = 1),
    data = dataLong
  )
coef <-
  cbind(coef(multLogitVGM), rep(1:N_CR, length(coef(multLogitVGM)) / N_CR))

# Calculate hazards

for(k in 1:N_CR) {
  coef1 <- coef[coef[, 2] == k, 1]
  temp <-
    coef1[[14]] * (data$GENDER == '2:Female') +
    coef1[[15]] * (data$RACE2 == '1') +
    coef1[[16]] * (data$RACE2 == '2') +
    coef1[[17]] * (data$ctrt == '2:RV to PA') +
    coef1[[18]] * (data$BWT) +
    coef1[[19]] * (data$PRENATDX == '1:Yes') +
    coef1[[20]] * (data$ATRESIA == '1:Yes') +
    coef1[[21]] * (data$OBSTRUCT == '1:Yes') +
    coef1[[22]] * (data$nordis_age)
  res <- matrix(NA, ncol = 13, nrow = nrow(data))
  res_prob <- matrix(NA, ncol = 13, nrow = nrow(data))
  for (i in 1:13) {
    if (i == 1) {
      res[, i] <- coef1[[i]] + temp
    }
    else {
      res[, i] <- coef1[[i]] + coef1[[1]] + temp
    }
  }
  for (i in 1:nrow(data)) {
    for (j in 1:13) {
      res_prob[i, j] <- exp(res[i, j])
    }
  }
 
  if(k==1) {
    c1 <- res_prob
  }  # c1 for S2P:Non-elective
  if (k == 2) {
    c2 <- res_prob
  }  # c2 for S2P:elective
  if (k == 3) {
    c3 <- res_prob
  }  # c3 for Death+Hrt Trans
}



h1 <- c1 * 1/(1 + c1 + c2 + c3) 
h2 <- c2 * 1/(1 + c1 + c2 + c3) 
h3 <- c3 * 1/(1 + c1 + c2 + c3) 

h0 <-  1/(c1 + c2 + c3 + 1)

# S(t): survival function 

St <- cumprod(h0[1,])   

# S(t-1)

St_1 <- c(1, St[1:12])   

# Generlized propensity score:

f_t1 <- matrix(NA, nrow = nrow(h1), ncol = 13)
for (i in 1:nrow(h1)) {
  f_t1[i, ] <- h1[i, ] * St_1
}

f_t2 <- matrix(NA, nrow = nrow(h2), ncol = 13)
for (i in 1:nrow(h2)) {
  f_t2[i, ] <- h2[i, ] * St_1
}

f_t3 <- matrix(NA, nrow = nrow(h3), ncol = 13)
for (i in 1:nrow(h3)) {
  f_t3[i, ] <- h3[i, ] * St_1
}

apply(f_t1, 2, mean)
apply(f_t2, 2, mean)
apply(f_t3, 2, mean)

avg_ft           <- cbind(apply(f_t1, 2, mean), apply(f_t2, 2, mean), apply(f_t3, 2, mean))
colnames(avg_ft) <- c('S2P:Non-Selective', 'S2P:Selective', 'Death_HrtTrans')
rownames(avg_ft) <- c(1:12, '>12')
avg_ft           <- avg_ft / sum(avg_ft)


# Plot FigureA3.

par(mar = c(5.1, 5.1, 4.1, 2.1))

plot(
  1:6,
  c(sum(avg_ft[1:3, 1]), avg_ft[4:7, 1], sum(avg_ft[8:13, 1])) ,
  type = 'l',
  lty = 4 ,
  col = 'blue',
  ylab = 'Probability',
  xlab = 'Norwood-S2P Interval Z (month)',
  cex.lab = 2.2,
  cex.axis = 2.1,
  ylim = c(0, 0.4),
  lwd = 2.8,
  xaxt = 'n'
)

axis(
  1,
  at = 1:6,
  labels = c('1-3', '4', '5', '6', '7', '>=8'),
  cex.axis = 2.1
)

lines(
  1:6,
  c(sum(avg_ft[1:3, 2]), avg_ft[4:7, 2], sum(avg_ft[8:13, 2])) ,
  type = 'l',
  col = 'red',
  lty = 1,
  lwd = 2.8
)

lines(
  1:6,
  c(sum(avg_ft[1:3, 3]), avg_ft[4:7, 3], sum(avg_ft[8:13, 3])) ,
  type = 'l',
  col = 'black',
  lty = 3,
  lwd = 2.8
)

legend(
  2.8,
  0.35,
  c(
    'R=1: Elective S2P',
    'R=2: Non-elective S2P',
    'R=3: Death/Heart Transplant'
  ),
  col = c('red', 'blue', 'black'),
  lty = c(1, 4, 3),
  cex = 1.6,
  lwd = 2.8
)

















