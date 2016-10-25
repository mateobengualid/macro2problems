# Cake eating, infinite time. Unless given, assume normalized cake.
# Since problem starts at t=0
cake_eating <- function(beta, x0=1.0) {
  # States. In R, everything starts at t=1. Extra niussance. -_-
  # Assume "infinite" equals T=100
  T = 100
  result <- data.frame(t=0:T, X=rep(x0,T+1), C=rep(0,T+1))
  for (i in 1:T) {
    result$X[i+1] <- result$X[i] * beta
  }
  
  # Controls.
  result$C <- result$X * (1-beta)
  
  return(result)
}

plot_ce <- function(ce, plot_title) {
  plot(x=ce$t, y=ce$X,type="l",ylim=c(0.0,1.05),xlab="t",ylab="Y",
       xaxs="i",yaxs="i",xaxt="n",lwd=2,main=plot_title)
  lines(x=ce$t, y=ce$C,col="red",lwd=2)
  axis(1, at = ce$t)
  
  # Legends have HUGE plot boxes. Plot the box by hand.
  leg <- function(pl) {return(legend(x="topright",legend=c("X", "C"),lwd=c(2.5,2.5),
                                     col=c("black", "red"),lty=c(1,1),cex=0.5,
                                     y.intersp=0.4,plot=pl, bty = "n"))}
  a <- leg(F)
  a <- a$rect
  mid <- a$top - 0.5*a$h
  reduction <- 0.5
  rect(xleft=a$left, ytop=mid+0.5*reduction*a$h, xright=a$left+a$w, ybottom=mid-0.5*reduction*a$h)
  leg(T)
}

# Plot.
par(mfrow=c(2,2))

# First example, small discount.
plot_ce(cake_eating(beta = 0.99), "B=0.99")

# Second example, low discount.
plot_ce(cake_eating(beta = 0.95), "B=0.95")

# Third example, mid discount.
plot_ce(cake_eating(beta = 0.85), "B=0.85")

# Fourth example, high discount.
plot_ce(cake_eating(beta = 0.5), "B=0.5")