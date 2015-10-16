library(package = "mFilter")

alpha <- 0.36
beta <- 0.99
gamma <- 1.7
delta <- 0.0012
ro <- 0.95
lambda <- delta ** (-delta)
sigma <- 0.01

# L is a constant.
L <- (1-alpha)*(beta*delta-beta+1)/
  ((1-beta+beta*delta-alpha*beta*delta)*gamma+(1-alpha)*(beta*delta-beta+1))
Anss <- 1.0
Knss <- (lambda*(alpha*beta*delta/(1-beta+beta*delta))**delta*(L**(1-alpha))**delta)**(1.0/(alpha-alpha*delta)) 

iterate <- function(df) {
  for(i in seq(nrow(df))[-201]) {
    # Stochastic version.
    df$A[i+1] <- df$A[i]**ro * exp(df$e[i])
    df$log_A[i+1] <- log(ro*df$A[i+1])
    df$Y[i] <- df$A[i] * (df$K[i]**alpha) * L**(1-alpha)
    df$C[i] <- df$Y[i] * (1-beta+beta*delta-alpha*beta*delta) / (beta*delta-beta+1)
    df$I[i] <- df$Y[i] - df$C[i]
    df$K[i+1] <- lambda * df$K[i]**(1-delta)*df$I[i]**delta
    
    # Nonstochastic version.
    df$Y_ns[i] <- (df$K_ns[i]**alpha) * L**(1-alpha)
    df$C_ns[i] <- df$Y_ns[i] * (1-beta+beta*delta-alpha*beta*delta) / (beta*delta-beta+1)
    df$I_ns[i] <- df$Y_ns[i] - df$C_ns[i]
    df$K_ns[i+1] <- lambda * df$K_ns[i]**(1-delta)*df$I_ns[i]**delta
  }
  i <- 201
  df$Y[i] <- df$A[i] * (df$K[i]**alpha) * L**(1-alpha)
  df$C[i] <- df$Y[i] * (1-beta+beta*delta-alpha*beta*delta) / (beta*delta-beta+1)
  df$I[i] <- df$Y[i] - df$C[i]
  
  df$Y_ns[i] <- (df$K_ns[i]**alpha) * L**(1-alpha)
  df$C_ns[i] <- df$Y_ns[i] * (1-beta+beta*delta-alpha*beta*delta) / (beta*delta-beta+1)
  df$I_ns[i] <- df$Y_ns[i] - df$C_ns[i]
  
  return(df)
}

# This list of gaussian values with mean 0 and sd 0.02 was generated from:
#https://www.random.org/gaussian-distributions/?num=201&mean=0.0&stdev=0.01&dec=10&col=1&notation=scientific&format=plain&rnd=id.mnb
exercise1f <- normal_0_001 <- read.csv("~/Projects/facultad/MACRO2/TP3/normal_0_001.csv")
exercise1f$A <- rep.int(Anss, nrow(exercise1f))
exercise1f$log_A <- rep.int(log(Anss), nrow(exercise1f))
exercise1f$K <- rep.int(Knss, nrow(exercise1f))
exercise1f$K_ns <- rep.int(Knss, nrow(exercise1f))
exercise1f <- iterate(exercise1f)
exercise1f$log_Y <- log(exercise1f$Y)
exercise1f$log_Y_ns <- log(exercise1f$Y_ns)
hp <- hpfilter(x=exercise1f$log_Y, freq = 1600, type = "lambda",drift = FALSE)

# Plot everything.
par(mfrow=c(1,1))
# Trend
# plot(exercise1f$t_plus_1,exercise1f$log_Y,type="l",xlab="t",ylab="log Y",xaxs="i",yaxs="i",lwd=2)
# lines(exercise1f$t_plus_1,exercise1f$log_Y_ns,col="red",lwd=2)
# lines(exercise1f$t_plus_1,hp$trend,col="green",lwd=2)
# legend(x="bottomright",legend=c("logY", "SS","Trend"),lwd=c(2.5,2.5,2.5),col=c("black", "red","green"), lty=c(1,1,1))
# Deviations and cyclic component
plot(exercise1f$t_plus_1,exercise1f$log_Y - exercise1f$log_Y_ns,type="l",xlab="t",ylab="log Y",
    ylim=c(-0.1,0.1), xaxs="i",yaxs="i",lwd=2)
abline(h=0,col="red",lwd=2)
lines(exercise1f$t_plus_1,hp$cycle,col="green",lwd=2)
legend(x="bottomleft",legend=c("Deviation from SS","Stationary State","Cyclic component"),lwd=c(2.5,2.5,2.5),col=c("black","red","green"), lty=c(1,1,1),cex = 0.8)