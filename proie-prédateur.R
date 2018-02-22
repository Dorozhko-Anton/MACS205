a = 0.25
b = 0.1
c = 1
d = 0.1

library(deSolve)

model <- function(time, y, parms) {
  dy1 <- y[1]*(a-b*y[2])
  dy2 <- -y[2]*(c-d*y[1])
  list(c(dy1, dy2))
}

yini = c(y1 = 20, y2 = 30)
times = seq(0, 100, 0.01)


out <- ode(times=times, y=yini, func=model, parms=NULL)

plot(out)

library(scatterplot3d)
scatterplot3d(out)


A = matrix(  c(0, 0, -1, 1, 0, -1, 0, 1, -1), ncol=3)

c(0, 0, cos(t))

expm1(A)
integrate()