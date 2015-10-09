exercise3d <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
alpha <- 0.5
g <- 0.02
sigma <- 0.02
beta <- 0.95
A0 <- 10
L <- 1.0
Knss <- A0 * L * (alpha * beta / exp(g))**(1.0/(1.0-alpha))

exercise3d$white_noise_mean_0_sd_d02 <- rep.int(0.0, nrow(exercise3d))
exercise3d$white_noise_mean_0_sd_d02[1] <- 2*sigma

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

exercise3d$ln_Y <- log(exercise3d$Y)
exercise3d$ln_y <- log(exercise3d$y)

# Plot everything.
par(mfrow=c(1,2))
plot(exercise3d$t[1:15],exercise3d$ln_y[1:15],type="l", main="3.d.I",ylim=c(-0.80,-0.72),xlab="t",ylab="ln y",xaxs="i",yaxs="i",lwd=2)
plot(exercise3d$t[1:15],exercise3d$ln_Y[1:15],type="l", main="3.d.II",xlab="t",ylab="ln Y",xaxs="i",yaxs="i",lwd=2)