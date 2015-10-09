alpha <- 0.5
beta <- 0.952
ro <- 0.5
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

exercise1c <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
exercise1c$tetha <- rep.int(tetha0, nrow(exercise1c))
exercise1c$log_tetha <- rep.int(log(tetha0), nrow(exercise1c))
exercise1c$k <- rep.int(k0, nrow(exercise1c))
exercise1c$k_ns <- rep.int(k0, nrow(exercise1c))
exercise1c <- iterate(exercise1c)

exercise1fcI <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
ro <- 0.1
exercise1fcI$tetha <- rep.int(tetha0, nrow(exercise1fcI))
exercise1fcI$log_tetha <- rep.int(log(tetha0), nrow(exercise1fcI))
exercise1fcI$k <- rep.int(k0, nrow(exercise1fcI))
exercise1fcI$k_ns <- rep.int(k0, nrow(exercise1fcI))
exercise1fcI <- iterate(exercise1fcI)

exercise1fcII <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
ro <- 0.9
exercise1fcII$tetha <- rep.int(tetha0, nrow(exercise1fcII))
exercise1fcII$log_tetha <- rep.int(log(tetha0), nrow(exercise1fcII))
exercise1fcII$k <- rep.int(k0, nrow(exercise1fcII))
exercise1fcII$k_ns <- rep.int(k0, nrow(exercise1fcII))
exercise1fcII <- iterate(exercise1fcII)

# Plot everything.
par(mfrow=c(1,3))
plot(exercise1c$t,exercise1c$y,type="l",xlab="t",ylab="y",
     main="1.c ro=0.5",ylim=c(0.3, 0.7),xaxs="i",yaxs="i",lwd=2)
lines(exercise1c$t,exercise1c$y_ns,col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))
plot(exercise1fcI$t,exercise1fcI$y,type="l",xlab="t",ylab="y",
     main="1.c.I ro=0.1",ylim=c(0.3, 0.7),xaxs="i",yaxs="i",lwd=2)
lines(exercise1fcI$t,exercise1fcI$y_ns,col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))
plot(exercise1fcII$t,exercise1fcII$y,type="l",xlab="t",ylab="y",
     main="1.c.II ro=0.9",ylim=c(0.3, 0.7),xaxs="i",yaxs="i",lwd=2)
lines(exercise1fcII$t,exercise1fcII$y_ns,col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))