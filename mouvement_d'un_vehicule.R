
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

plot(x, y, type='l')

#out <- ode(times=times, y=yini, func=model, parms=NULL)


runge_kutta = function(times, y, func, parms) {
  state = y
  
  trace = list(c(state))
  
  for (i in seq(1, length(times)-1, 1)) {
    h = times[i+1] - times[i]
    
    p1 = func(times[i], state, parms)[[1]]
    t2 = times[i] + 0.5*h
    y2 = state + 0.5*h*p1
    p2 = func(t2, y2, parms)[[1]]
    
    t3 = t2
    y3 = state + 0.5*h*p2
    p3 = func(t3, y3, parms)[[1]]
    
    t4 = times[i] + h
    y4 = state + h*p3
    p4 = func(t4, y4, parms)[[1]]
    
    state = state + h/6*(p1 + 2*p2 + 2*p3 + p4)
    
    trace = c(trace, list(c(state)))
  }
  trace
}

out_rk <- runge_kutta(times=times, y=yini, func=model, parms = NULL)

x = c()
y = c()
for (i in seq(1, length(out), 1)) {
  x = c(x, out_rk[[i]][1])
  y = c(y, out_rk[[i]][2])
}

lines(x, y, type='l', col='red')
