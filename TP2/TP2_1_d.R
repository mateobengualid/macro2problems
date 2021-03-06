exercise1c <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
alpha <- 0.5
beta <- 0.952
ro <- 0.5
sigma <- 0.02
k0 <- 0.005
tetha0 <- 1.0

exercise1c$tetha <- rep.int(tetha0, nrow(exercise1c))
exercise1c$log_tetha <- rep.int(log(tetha0), nrow(exercise1c))
exercise1c$k <- rep.int(k0, nrow(exercise1c))
exercise1c$k_ns <- rep.int(k0, nrow(exercise1c))

# The only "hard" to write are self-referenced variables.
for(i in seq(nrow(exercise1c))[-101]) {
  # Stochastic version.
  exercise1c$log_tetha[i+1] <- ro*exercise1c$log_tetha[i] + exercise1c$white_noise_mean_0_sd_d02[i]
  exercise1c$tetha[i+1] <- exp(exercise1c$log_tetha[i+1])
  exercise1c$y[i] <- exercise1c$tetha[i] * (exercise1c$k[i]**alpha)
  exercise1c$c[i] <- exercise1c$y[i] * (1-alpha*beta)  
  exercise1c$k[i+1] <- exercise1c$y[i] - exercise1c$c[i]
  
  # Nonstochastic version.
  exercise1c$y_ns[i] <- exercise1c$k_ns[i]**alpha
  exercise1c$c_ns[i] <- exercise1c$y_ns[i] * (1-alpha*beta)
  exercise1c$k_ns[i+1] <- exercise1c$y_ns[i] - exercise1c$c_ns[i]
}

# Plot everything.
par(mfrow=c(1,1))
plot(exercise1c$t[1:51],exercise1c$y[1:51],type="l",
     main="1.d",xlab="t",ylab="y",ylim=c(0.0, 0.52),xaxs="i",yaxs="i",lwd=2)
lines(x = exercise1c$t[1:51], y = exercise1c$y_ns[1:51],col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))