exercise1c <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
alpha <- 0.5
beta <- 0.952
ro <- 0.5
sigma <- 0.02
k0 <- (alpha*beta)**(1.0/(1.0-alpha))
tetha0 <- 1.0

# This time, €1 (the first one) is 2*sigma and 0 afterwards.
exercise1c$white_noise_mean_0_sd_d02 <- rep.int(0.0, nrow(exercise1c))
exercise1c$white_noise_mean_0_sd_d02[1] <- 2 * sigma
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

# Non self-referenced variables are easier to write.
exercise1c$y <- 
 


# Plot everything.
plot(exercise1c$t[1:26],exercise1c$y[1:26],type="l",xlab="t",ylab="y",
     main="1.e",ylim=c(0.4, 0.52),xaxs="i",yaxs="i",lwd=2)
lines(exercise1c$t[1:26],exercise1c$y_ns[1:26],col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))