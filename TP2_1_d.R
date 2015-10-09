exercise1c <- read.csv("~/Projects/facultad/MACRO2/normal_0_002.csv")
alpha <- 0.5
beta <- 0.952
ro <- 0.5
sigma <- 0.02
k0 <- 0.005
tetha0 <- 1.0

exercise1c$tetha <- rep.int(tetha0, nrow(exercise1c))
exercise1c$k <- rep.int(k0, nrow(exercise1c))
exercise1c$k_ns <- rep.int(k0, nrow(exercise1c))

# The only "hard" to write are self-referenced variables.
for(i in seq(nrow(exercise1c))[-101]) {
  exercise1c$tetha[i+1] <-
    (exercise1c$tetha[i] ** ro) * exp(exercise1c$white_noise_mean_0_sd_d02[i])
  
  exercise1c[i+1,"k"] <- alpha * beta * exercise1c[i,"tetha"] * exercise1c[i+1,"k"] ** alpha
  exercise1c[i+1,"k_ns"] <- alpha * beta * exercise1c[i+1,"k_ns"] ** alpha
}

# Non self-referenced variables are easier to write.
exercise1c$y <- exercise1c$tetha * (exercise1c$k**alpha)
exercise1c$c <- exercise1c$y * (1-alpha*beta)
exercise1c$y_ns <- exercise1c$k_ns**alpha
exercise1c$c_ns <- exercise1c$y_ns * (1-alpha*beta)

# Plot everything.
plot(exercise1c$t[0:50],exercise1c$y[0:50],type="l",xlab="t",ylab="y",ylim=c(0.0, 0.52),xaxs="i",yaxs="i",lwd=2)
lines(exercise1c$t[0:50],exercise1c$y_ns[0:50],col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))