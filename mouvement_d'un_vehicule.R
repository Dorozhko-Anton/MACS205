
x0 = 0
y0 = 0
v0 = 0
phi0 = 0


a = function(t) {
  if (t > 0 && t < 3.4) {
    if (t > 1 && t < 2.7) {
      return (-1)
    } else {
      return (1)
    }
  }
  return (0)
}

phi = function(t) {
  if (t > 0 && t < 1) return (1)
  if (t > 1 && t < 2) return (-1)
  return (0)
}

library(deSolve)

model = function(time, state, parms) {
  x = state[1]
  y = state[2]
  v = state[3]
  theta = state[4]
  
  dx <- v*cos(theta)
  dy <- v*sin(theta)
  dv <- a(time)
  dtheta <- phi(time)
  list(c(dx, dy, dv, dtheta))
}

euler = function(times, y, func, parms) {
  state = y
  
  trace = list(c(state))
  
  for (i in seq(1, length(times)-1, 1)) {
    h = times[i+1] - times[i]
    print(func(times[i], state, parms)[[1]])
    state = state + func(times[i], state, parms)[[1]] * h
    trace = c(trace, list(c(state)))
  }
  trace
}


yini = c(x = 0, y = 0, v=0, phi=0)
times = seq(0, 10, 0.01)

out <- euler(times=times, y=yini, func=model, parms = NULL)

x = c()
y = c()
for (i in seq(1, length(out), 1)) {
  x = c(x, out[[i]][1])
  y = c(y, out[[i]][2])
}

plot(x, y)

#out <- ode(times=times, y=yini, func=model, parms=NULL)



