alpha <- 0.36
beta <- 0.99
gamma <- 1.7
delta <- 0.0012
ro <- 0.95
lambda <- delta ** (-delta)

# L is a constant.
L <- (1-alpha)*(beta*delta-beta+1)/
  ((1-beta+beta*delta-alpha*beta*delta)*gamma+(1-alpha)*(beta*delta-beta+1))
Anss <- 1.0
Knss <- (lambda*(alpha*beta*gamma/(1-beta+beta*delta))**delta*(L**(1-alpha))**delta)**(1.0/(alpha-alpha*delta)) 

iterate <- function(df) {
  for(i in seq(nrow(df))[-201]) {
    # Stochastic version.
    df$A[i+1] <- df$A[i]**ro * exp(df$e[i])
    df$log_A[i+1] <- log(ro*df$A[i+1])
    df$Y[i] <- df$A[i] * (df$K[i]**alpha) * L**(1-alpha)
    df$C[i] <- df$Y[i] * (beta*delta-beta+1) / (1-beta+beta*delta-alpha*beta*delta)
    df$I[i] <- df$Y[i] - df$C[i]
    df$K[i+1] <- lambda * df$K[i]**(1-delta)*df$I[i]**delta
    
    # Nonstochastic version.
    df$Y_ns[i] <- (df$K_ns[i]**alpha) * L**(1-alpha)
    df$C_ns[i] <- df$Y_ns[i] * (beta*delta-beta+1) / (1-beta+beta*delta-alpha*beta*delta)
    df$I_ns[i] <- df$Y_ns[i] - df$C_ns[i]
    df$K_ns[i+1] <- lambda * df$K_ns[i]**(1-delta)*df$I_ns[i]**delta
  }
  
  return(df)
}

# This list of gaussian values with mean 0 and sd 0.02 was generated from:
#https://www.random.org/gaussian-distributions/?num=201&mean=0.0&stdev=0.01&dec=10&col=1&notation=scientific&format=plain&rnd=id.mnb
exercise1e <- normal_0_001 <- read.csv("~/Projects/facultad/MACRO2/TP3/normal_0_001.csv")
exercise1e$A <- rep.int(Anss, nrow(exercise1e))
exercise1e$log_A <- rep.int(log(Anss), nrow(exercise1e))
exercise1e$K <- rep.int(Knss, nrow(exercise1e))
exercise1e$K_ns <- rep.int(Knss, nrow(exercise1e))
exercise1e <- iterate(exercise1e)

# Plot everything.
par(mfrow=c(1,1))
plot(exercise1e$t_plus_1,exercise1e$Y,type="l",xlab="t",ylab="y",
     main="1.c ro=0.5",xaxs="i",yaxs="i",lwd=2)
lines(exercise1e$t_plus_1,exercise1e$Y_ns,col="red",lwd=2)
legend(x="bottomright",legend=c("Shocks", "No Shocks"),lwd=c(2.5,2.5),col=c("black", "red"), lty=c(1,1))