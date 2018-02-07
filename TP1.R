# Ex1

# Ex2

x0 = log(6/5)

Gn = function(x0, n) {
  x <- c(x0)

  for (i in (2:n)) {
    x[i] <- 1/(i-1) - x[i-1]          
  }
  return (x)
}

# Ex3

epsilon = function() {
  e <- 1
  while ((1+e > 1)) {
    e <- e / 2
  }
   
  return (2*e)  
}

# Ex4

# approximation 1
myDiff1 = function(x) {
  if (x == 0) return (1)
  else return ( (exp(x) - 1)/x)
}

# approximation 2
myDiff2 = function(x) {
  y = exp(x)
  if (y==1) return (1)
  else return ( (y-1)/log(y) )
}

plot_instability = function() {
  y1 <- c()
  y2 <- c()
  eps = .Machine$double.eps/2
  
  x = matrix(c(-10:10)*eps)
  
  y1 <- apply(x, 1, myDiff1)
  y2 <- apply(x, 1, myDiff2)
  

  plot(x, y1, type="l")
  lines(x, y2, col="green")
  abline(v=.Machine$double.eps, col="red")
}
