# HP Filter library.
library(package = "mFilter")

# Calibration.
theta <- 0.36 #Kidland & Prescott.
delta <- 0.025 # 10% anual. Good compromise for different rates of depr.
beta <- 0.99 # 4% real interest rate.
A <- 2 # This results in Hnss of 0.33 in steady state. So, 8 hours a day.
H0 <- 0.53 # Result from previous values.
sigma_e <- 0.00712 # From GDP for the U.S. Economy.
gamma <- 0.95 # "A more detailed study of the statistical properties of
# this technology shock is planned but has not yet been carried out."
totalPeriods <- 115 # Postwar quarters.
totalSimulations <- 100 # As per Hansen's paper.

# NSS variables.
lambdaNSS <- 1 # Unconditional mean.
rNSS <- (1 - beta + delta*beta)/beta
wNSS <- (1-theta)*lambdaNSS*(rNSS/theta/lambdaNSS)**(-theta/(1-theta))
KNSS <- theta*wNSS /((A+1-theta)*rNSS -A*delta*theta)
CNSS <- (rNSS/theta - delta)*KNSS
HNSS <- KNSS * (rNSS/lambdaNSS/theta)**(1/(1-theta))
INSS <- delta * KNSS
YNSS <- rNSS*KNSS/theta

# Xt = [Kt+1]
# (Yt)' <- [Yt Ct It Ht rt wt]
# Zt <- [LAMBDAt]
# et+1 <- [et+1]

# Then:
# 0 = A Xt + B Xt-1 + C Yt + D Zt
# 0 = Et {F Xt+1 + G Xt + H Xt-1 + J yt+1 + K Yt, L Zt+1 + M Zt}
# zt+1 = N zt + et+1, E[e t+1] = 0

AA = matrix(c(0, 0, 0, 0, 1, 0))

BB = matrix(c(-1, theta, 0, 0, -(1-delta), 0))

CC = rbind(c( 1,      -1, 1, -YNSS,      0,              0), # Yt)
           c( 0,       0, 0, +CNSS,      0,             -1), # Ct)
           c( 0,       0, 0, +INSS, -delta,              0), # It)
           c( 0, 1-theta,-1,     0,      0, -HNSS/(1-HNSS)), # Ht)
           c(-1,       0, 0,     0,      0,              0), # rt)
           c( 0,       0,-1,     0,      0,              1)) # wt)
CC <- t(CC) # Get Equations in the rows.

DD = matrix(c( 0, 1, 0, 0, 0, 0))

FF = matrix(c(0))

GG = matrix(c( 0 ))

HH = matrix(c( 0 ))

JJ = t(matrix(c(0, -1, 0, 0, beta*rNSS, 0)))

KK = t(matrix(c(0, 1, 0, 0, 0, 0)))

LL = matrix(c( 0 ))

MM = matrix(c( 0 ))

NN = matrix(c(gamma))

Sigma = matrix(c( sigma_e^2  ))

# Setting the options:
l_equ <- dim(AA)[1]
m_states <-dim(AA)[2]
n_endog <-dim(CC)[2]
k_exog <-dim(DD)[2]
VARNAMES <- matrix(c('K', 'Y','C','I','H','r','w', 'lambda'))
PERIOD     <- 4 # number of periods per year, i.e. 12 for monthly, 4 for quarterly
GNP_INDEX  <- 4 # Index of output among the variables selected for HP filter
IMP_SELECT <- 1:(m_states+n_endog+k_exog)
HP_SELECT  <- 1:(m_states+n_endog+k_exog) # Selecting the variables for the HP Filter calcs.
DO_SIMUL   <- 0 # Calculates Simulations
DO_MOMENTS <- 0 # Calculates Moments
DO_PLOTS <- 0 # Plots
DISPLAY_AT_THE_END <- 0 # Stop filling my console with text!
DISPLAY_IMMEDIATELY <- 0 # Stop asking for bloody enters!
DO_QZ <- 0 # Solve using generalized eigenvalues instead of Schur's factorization.
source("./UhligR/do_it.R", chdir = TRUE)

# Iterate over state and control variables. 
iterate <- function(shocks, k0 = 0, lambda0 = 0, periods = totalPeriods) {
  states <- data.frame(k=rep_len(k0, periods),
                       lambda=rep_len(lambda0, periods),
                       t=0:(totalPeriods-1),
                       e1=shocks)
  
  # Iterate for states.
  for(i in 1:(periods-1)) {
    # X(t) = PP %*% X(t-1) + QQ %*% Z(t).
    states$k[i+1] <- PP * states$k[i] + QQ * states$lambda[i]
    states$lambda[i+1] <- gamma * states$lambda[i] + shocks[i]
  }
  
  # Solve for controls.
  # Y(t) = RR %*% X(t-1) + SS %*% Z(t)
  controls <- t(RR %*% states$k + SS %*% states$lambda)
  colnames(controls)<-c('y','c','i','h','r','w')
  controls <- data.frame(controls)
  controls$t <- states$t
    
  return(merge.data.frame(states,controls, by = "t"))
}

# Iterate.
stationary <- iterate(shocks = rep.int(0,totalPeriods))
simulation <- iterate(shocks = rnorm(totalPeriods, mean = 0, sd=sigma_e*100))

# Recover the log(X) from the difference between logarithm and log(Xnss)
simulation$log_K <- simulation$k + log(KNSS)
simulation$log_LAMBDA <- simulation$lambda + log(lambdaNSS)
simulation$log_Y <- simulation$k + log(YNSS)
simulation$log_C <- simulation$c + log(CNSS)
simulation$log_I <- simulation$i + log(INSS)
simulation$log_H <- simulation$h + log(HNSS)
simulation$log_R <- simulation$r + log(rNSS)
simulation$log_W <- simulation$w + log(wNSS)

hpcycle <- function(x) hpfilter(x=x, freq = 1600, type = "lambda",drift = FALSE)$cycle

par(mfrow=c(2,4))
plot(stationary$t,hpcycle(simulation$k),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Capital (K)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,hpcycle(simulation$lambda),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Technology (Lambda)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,hpcycle(simulation$y),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Product (Y)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,hpcycle(simulation$c),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Consumption (C)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,hpcycle(simulation$i),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Investment (I)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,hpcycle(simulation$h),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Labor (H)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,hpcycle(simulation$r),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Interest rate (r)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,hpcycle(simulation$w),col="red",type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,main="Wage (w)")
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")

# Calculate, for Y,C,I,K,H and w
# a) mean percentage deviation from trend (with standard error)
# b) correlation against output/production (with standard error)
table <- data.frame(mpdft_Y = rep(0, totalSimulations), cor_Y = rep(0, totalSimulations),
                    mpdft_C = rep(0, totalSimulations), cor_C = rep(0, totalSimulations),
                    mpdft_I = rep(0, totalSimulations), cor_I = rep(0, totalSimulations),
                    mpdft_K = rep(0, totalSimulations), cor_K = rep(0, totalSimulations),
                    mpdft_H = rep(0, totalSimulations), cor_H = rep(0, totalSimulations),
                    mpdft_w = rep(0, totalSimulations), cor_w = rep(0, totalSimulations))
for(i in 1:totalSimulations) {
  # Simulate once.
  sh <- rnorm(totalPeriods, mean = 0, sd=sigma_e*100)
  sim <- iterate(shocks = sh)
  
  # Detrend.
  sim$cycle_Y <- hpcycle(sim$y)
  sim$cycle_C <- hpcycle(sim$c)
  sim$cycle_I <- hpcycle(sim$i)
  sim$cycle_K <- hpcycle(sim$k)
  sim$cycle_H <- hpcycle(sim$h)
  sim$cycle_w <- hpcycle(sim$w)
  
  # Get the percentage deviation from trend
  table$mpdft_Y[i] <- sd(sim$cycle_Y)
  table$mpdft_C[i] <- sd(sim$cycle_C)
  table$mpdft_I[i] <- sd(sim$cycle_I)
  table$mpdft_K[i] <- sd(sim$cycle_K)
  table$mpdft_H[i] <- sd(sim$cycle_H)
  table$mpdft_w[i] <- sd(sim$cycle_w)
  
  # Get the correlations.
  table$cor_Y[i] <- cor(x = sim$cycle_Y, y = sim$cycle_Y)
  table$cor_C[i] <- cor(x = sim$cycle_Y, y = sim$cycle_C)
  table$cor_I[i] <- cor(x = sim$cycle_Y, y = sim$cycle_I)
  table$cor_K[i] <- cor(x = sim$cycle_Y, y = sim$cycle_K)
  table$cor_H[i] <- cor(x = sim$cycle_Y, y = sim$cycle_H)
  table$cor_w[i] <- cor(x = sim$cycle_Y, y = sim$cycle_w)
  
  if(i %% 10 == 0) {
    print(i)
  }
}
print("Divisible labor")
for(v in colnames(table)) {
  print(paste(v,": ", mean(table[,v]), " (",sd(table[,v]),")"))
}
