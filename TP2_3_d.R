exercise3d <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
alpha <- 0.5
g <- 0.02
beta <- 0.95
A0 <- 10
L <- 1.0
Knss <- (alpha * beta / exp(g))**(1.0/(1.0-alpha))

exercise3d$A <- rep.int(A0, nrow(exercise3d))
exercise3d$K <- rep.int(Knss, nrow(exercise3d))

for(i in seq(nrow(exercise3d))[-101]) {
  # Big version.
  exercise3d$A_mult[i] <- exp(g+exercise3d$white_noise_mean_0_sd_d02[i])
  exercise3d$A[i+1] <- exercise3d$A[i] * exercise3d$A_mult[i]
  exercise3d$Y[i] <- exercise3d$K[i] ** alpha * (L * exercise3d$A[i]) ** (1 - alpha)
  exercise3d$C[i] <- (1 - alpha * beta) * L ** (1-alpha) * exercise3d$A[i] ** (1-alpha) * exercise3d$K[i] ** alpha
  exercise3d$K[i+1] <- exercise3d$Y[i] - exercise3d$C[i]
  
  # Effective version.
  exercise3d$y[i] <- exercise3d$Y[i] / exercise3d$A[i] / L
  exercise3d$c[i] <- exercise3d$C[i] / exercise3d$A[i] / L
  exercise3d$k[i] <- exercise3d$K[i] / exercise3d$A[i] / L
}

# Plot everything.
par(mfrow=c(1,1))
plot(exercise3d$t,log(exercise3d$y),type="l",
     main="3.d",xlab="t",ylab="y",xaxs="i",yaxs="i",lwd=2)
lines(exercise3d$t,log(exercise3d$Y),col="red",lwd=2)
legend(x="bottom",legend=c("ln y", "ln Y"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))