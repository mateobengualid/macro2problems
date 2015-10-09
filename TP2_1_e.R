alpha <- 0.5
beta <- 0.952
sigma <- 0.02
k0 <- (alpha*beta)**(1.0/(1.0-alpha))
tetha0 <- 1.0

iterate <- function(df) {
  for(i in seq(nrow(df))[-101]) {
    # Stochastic version.
    df$log_tetha[i+1] <- ro*df$log_tetha[i] + df$white_noise_mean_0_sd_d02[i]
    df$tetha[i+1] <- exp(df$log_tetha[i+1])
    df$y[i] <- df$tetha[i] * (df$k[i]**alpha)
    df$c[i] <- df$y[i] * (1-alpha*beta)  
    df$k[i+1] <- df$y[i] - df$c[i]
    
    # Nonstochastic version.
    df$y_ns[i] <- df$k_ns[i]**alpha
    df$c_ns[i] <- df$y_ns[i] * (1-alpha*beta)
    df$k_ns[i+1] <- df$y_ns[i] - df$c_ns[i]
  }
  
  return(df)
}

# This time, €1 (the first one) is 2*sigma and 0 afterwards.
exercise1c <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
ro <- 0.5
exercise1c$white_noise_mean_0_sd_d02 <- rep.int(0.0, nrow(exercise1c))
exercise1c$white_noise_mean_0_sd_d02[1] <- 2 * sigma
exercise1c$tetha <- rep.int(tetha0, nrow(exercise1c))
exercise1c$log_tetha <- rep.int(log(tetha0), nrow(exercise1c))
exercise1c$k <- rep.int(k0, nrow(exercise1c))
exercise1c$k_ns <- rep.int(k0, nrow(exercise1c))
exercise1c <- iterate(exercise1c)

exercise1c2 <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
ro <- 0.1
exercise1c2$white_noise_mean_0_sd_d02 <- rep.int(0.0, nrow(exercise1c2))
exercise1c2$white_noise_mean_0_sd_d02[1] <- 2 * sigma
exercise1c2$tetha <- rep.int(tetha0, nrow(exercise1c2))
exercise1c2$log_tetha <- rep.int(log(tetha0), nrow(exercise1c2))
exercise1c2$k <- rep.int(k0, nrow(exercise1c2))
exercise1c2$k_ns <- rep.int(k0, nrow(exercise1c2))
exercise1c2 <- iterate(exercise1c2)

exercise1c3 <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
ro <- 0.9
exercise1c3$white_noise_mean_0_sd_d02 <- rep.int(0.0, nrow(exercise1c3))
exercise1c3$white_noise_mean_0_sd_d02[1] <- 2 * sigma
exercise1c3$tetha <- rep.int(tetha0, nrow(exercise1c3))
exercise1c3$log_tetha <- rep.int(log(tetha0), nrow(exercise1c3))
exercise1c3$k <- rep.int(k0, nrow(exercise1c3))
exercise1c3$k_ns <- rep.int(k0, nrow(exercise1c3))
exercise1c3 <- iterate(exercise1c3)

# Plot everything.
par(mfrow=c(1,3))
plot(exercise1c$t[1:26],exercise1c$y[1:26],type="l",xlab="t",ylab="y",
     main="1.e",ylim=c(0.4, 0.52),xaxs="i",yaxs="i",lwd=2)
lines(exercise1c$t[1:26],exercise1c$y_ns[1:26],col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))
plot(exercise1c2$t[1:26],exercise1c2$y[1:26],type="l",xlab="t",ylab="y",
     main="1.e.I ro=0.1",ylim=c(0.4, 0.52),xaxs="i",yaxs="i",lwd=2)
lines(exercise1c2$t[1:26],exercise1c2$y_ns[1:26],col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))
plot(exercise1c3$t[1:26],exercise1c3$y[1:26],type="l",xlab="t",ylab="y",
     main="1.e.II ro=0.9",ylim=c(0.4, 0.52),xaxs="i",yaxs="i",lwd=2)
lines(exercise1c3$t[1:26],exercise1c3$y_ns[1:26],col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))