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

# We introduced a B variable for ease. Calculate it now.
B <- -A*log(1-H0)/H0

# NSS variables.
lambdaNSS <- 1 # Unconditional mean.
rNSS <- (1 - beta + delta*beta)/beta
wNSS <- (1-theta)*lambdaNSS*(rNSS/theta/lambdaNSS)**(-theta/(1-theta))
KNSS <- theta*wNSS /(B*(rNSS-delta*theta))
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

CC = rbind(c( 1,      -1, 1, -YNSS,      0, 0), # Yt)
           c( 0,       0, 0, +CNSS,      0,-1), # Ct)
           c( 0,       0, 0, +INSS, -delta, 0), # It)
           c( 0, 1-theta,-1,     0,      0, 0), # Ht)
           c(-1,       0, 0,     0,      0, 0), # rt)
           c( 0,       0,-1,     0,      0, 1)) # wt)
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
impulse <- iterate(shocks = c(1, rep.int(0,totalPeriods-1)))
par(mfrow=c(2,4))
plot(stationary$t,stationary$k,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.1,1.2),main="Capital (K)")
lines(impulse$t,impulse$k,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,stationary$lambda,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.1,1.2),main="Technology (Lambda)")
lines(impulse$t,impulse$lambda,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,stationary$y,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.1,2),main="Product (Y)")
lines(impulse$t,impulse$y,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,stationary$c,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.1,1.0),main="Consumption (C)")
lines(impulse$t,impulse$c,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,stationary$i,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.25,7),main="Investment (I)")
lines(impulse$t,impulse$i,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,stationary$h,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.2,1.6),main="Labor (H)")
lines(impulse$t,impulse$h,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,stationary$r,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.5,2),main="Interest rate (r)")
lines(impulse$t,impulse$r,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")
plot(stationary$t,stationary$w,type="l",xlab="Qs",ylab="% deviation from SS",xaxs="i",yaxs="i",lwd=2,ylim=c(-0.1,1),main="Wage (w)")
lines(impulse$t,impulse$w,col="red",lwd=2)
legend(x="topright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1),cex=0.4,bty = "n")