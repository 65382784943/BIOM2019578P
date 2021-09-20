



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



# Sort by blind_id
# For each R, use right f_t1, f_t2 or f_t3, etc. 

f_t2           <- as.data.frame(f_t2)
f_t2           <- cbind(data$blind_id, f_t2)
colnames(f_t2) <- c('blind_id', 1:12, '>12')

f_t2           <- f_t2[order(f_t2$blind_id),]
data2          <- data2[order(data2$blind_id),]



# Get GPS for S2P patients

temp           <-
  right_join(f_t2, select(data2, blind_id), by = 'blind_id')
f_t2           <- temp
rm(temp)

prob_for6group            <-
  cbind(f_t2[, 1], apply(f_t2[, 2:4], 1, sum), (f_t2[, 5:8]), apply(f_t2[, 9:14], 1, sum))
colnames(prob_for6group)  <- c('blind_id', 1:6)
prob_for6group2           <-
  cbind(prob_for6group[, 1], sweep(prob_for6group[,-1], 1, rowSums(prob_for6group[,-1]), FUN =
                                     "/"))
colnames(prob_for6group2) <- c('blind_id', 1:6)
rm(prob_for6group)

STG2_month_grp            <-
  as.numeric(cut(
    data2$STG2_age - data2$norw_age,
    breaks = c(0, seq(91.5, 213.5, by = 30.5), 800),
    labels = c('1-3', 4:7, '>=8')
  ))

temp <-
  left_join(prob_for6group2,
            select(data2, blind_id, event2_time, censor2),
            by = 'blind_id')
temp[is.na(temp$event2_time), ]$event2_time <-
  ceiling(mean(temp[temp$censor2 == 0, ]$event2_time, na.rm = T))
temp <- cbind(temp, STG2_month_grp)
rm(STG2_month_grp, prob_for6group2)
temp <-
  data.frame(temp, wgt = 1 / apply(temp[, c(2:7, 10)], 1, function(x) {
    x[x[7]]
  }))

# Get stabilized weights

Marginal     <- apply(temp[, 2:7], 2, mean)
Stb_wgt      <- rep(0, nrow(temp))
for (i in 1:nrow(temp)) {
  Stb_wgt[i] <- Marginal[temp$STG2_month_grp[i]] * temp$wgt[i]
}
temp         <- data.frame(temp, Stb_wgt)
rm(Stb_wgt)

# Cox ph regression with IPSW

fit1 <-
  coxph(
    Surv(event2_time, censor2) ~ lsp(temp$STG2_month_grp, knot = c(3, 4)),
    weights = Stb_wgt,
    data = temp,
    x = T
  )
summary(fit1)

# Use AIC to select knots

extractAIC(fit1)

# Get covariance matrix for \alpha 's; also available by 'vcov(fit1)'

hr1       <- smoothHR(data = temp, coxfit = fit1)
temp_plot <-
  predict(
    hr1,
    predictor = "STG2_month_grp",
    prob = 0,
    conf.level = 0.95,
    prediction.values = 1:6,
    ref.label = "Z"
  )

# Plot Figure 3

Arrow     <-
  data.frame(avg = temp_plot[, 2],
             sdev = (temp_plot[, 4] - temp_plot[, 3]) / 2 / 1.96)
Arrow$low <- Arrow[, 1] - Arrow[, 2]
Arrow$up  <- Arrow[, 1] + Arrow[, 2]
x         <- 1:6


par(mar = c(5.1, 5.1, 4.1, 2.1))

attach(Arrow)

plot(
  x,
  avg,
  xlab = 'Norwood-S2P interval (month)',
  ylab = 'Ln Hazard Ratio (Z, Zref = 6)',
  cex.lab = 1.5,
  cex.axis = 1.5,
  lwd = 2,
  ylim = c(-0.2, 2),
  xaxt = 'n'
)

axis(
  1,
  at = c(1:6),
  labels = c('1-3', '4', '5', '6', '7', '>=8'),
  cex.axis = 1.6
)

arrows(
  x,
  low,
  x,
  up,
  length = 0.03,
  angle = 90,
  code = 3,
  lwd = 2
)

detach(Arrow)
abline(h = 0)














