# Cake eating, finite time. Unless given, assume normalized cake.
# Since problem starts at t=0, at T you have T+1 periods.
cake_eating <- function(beta, T, x0=1.0) {
  # States. In R, everything starts at t=1. Extra niussance. -_-
  result <- data.frame(t=0:T, X=rep(x0,T+1), C=rep(0,T+1))
  for (i in 1:T) {
    # Calculate the string of B+B^2+B^3+...BT-(i-1).
    # That (i-1) stems from R arrays start at 1, not at 0.
    beta_sum <- 0.0
    for (j in 1:(T-i+1)) {
      beta_sum <- beta_sum + beta**(j)
    }
    result$X[i+1] <- result$X[i] * beta_sum / (1+beta_sum)
    
    # It's not necessary to calculate the controls here,
    # but it's cheaper in computational time.
    result$C[i] <- result$X[i] / (1+beta_sum)
  }
  # The last cake eating ends up calculated outside the loop. Ah well...
  result$C[T+1] <- result$X[T+1]
  
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

# First example, 4 periods, no discount.
plot_ce(cake_eating(beta = 1.0, T = 3), "T=3, B=1.0")

# Second example, 4 periods, high discount.
plot_ce(cake_eating(beta = 0.50, T = 3), "T=3, B=0.5")

# Third example, 100 periods, no discount.
plot_ce(cake_eating(beta = 1.0, T = 99), "T=99, B=1.0")

# Fourth example, 100 periods, mild discount.
plot_ce(cake_eating(beta = 0.95, T = 99), "T=99, B=0.95")