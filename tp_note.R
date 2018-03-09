library(deSolve)

w0 = 1.2
eps = 0.2

eq_van_der_pol = function(time, state, parms) {
  x = state[1]
  z = state[2]
  
  dx <- z
  dz <- eps*w0*(1-x^2)*z - w0^2*x
  list(c(dx, dz))
}

euler = function(times, y, func, parms) {
  state = y
  
  trace = matrix(c(times[1], state), nrow=1)
  
  for (i in seq(1, length(times)-1, 1)) {
    h = times[i+1] - times[i]
    
    state = state + func(times[i], state, parms)[[1]] * h
    
    trace = rbind(trace, c(times[i+1], c(state)))
  }
  trace
}

yini = c(x = 1, z = 0)
times = seq(0, 50, 0.1)

out <- ode(times=times, y=yini, func=eq_van_der_pol, parms = NULL)
out_euler <- euler(times=times, y=yini, func=eq_van_der_pol, parms = NULL)

plot(out[,2], type='l')
lines(out_euler[,2], type='l', col='red')


plot(out[,3], type='l')
lines(out_euler[,3], type='l', col='red')
# 1.3

hs = c(1, 0.5, 0.1, 0.01, 0.001)

err_T = c()

for (h in hs) {
  yini = c(x=1, z=0)
  times = seq(0, 50, h)
  
  out <- ode(times=times, y=yini, func=eq_van_der_pol, parms = NULL)
  out_euler <- euler(times=times, y=yini, func=eq_van_der_pol, parms = NULL)
  
  err_T = c(err_T, abs(out[,2][50] - out_euler[,2][50]))
}

plot(log(hs), log(err_T), type='l')


# 1.4
yini = c(x=1, z=0)
times = seq(0, 50, 0.001)
system.time(out_euler <- euler(times=times, y=yini, func=eq_van_der_pol, parms = NULL))

F = function(x, z) {
  if (x >= 3 && x <= -3) {
    return(0)
  }
  if (z >= 3 && z <= -3) {
    return(0)
    }
  
  return (-abs(-2*x*eps*w0*z - w0^2 + eps*w0*(1-x^2)))
}

a=seq(-3, 3, by=0.001)
b=seq(-3, 3, by =0.001)
r = outer(a, b, FUN=F)


#2.2.1
a = 0.95
b = 0.4
tau = 11

Iext_2 = function(t) {
  (0.59)
}

eq_fitzh_nagumo= function(time, state, parms) {
  with(as.list(c(state, parms)), {
    dv <- v - v^3/3 - w + Iext(time)
    dw <- (v + a - b*w)/tau
    return(list(c(dv, dw)))
  })
}

parms <- c(
  a = a,
  b = b,
  tau = tau,
  Iext = Iext_2
)
yini = c(v=0, w=0)
times = seq(0, 300, 1)

out_2 <- euler(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)
plot(out_2[,2], type='l')
plot(out_2[,3], type='l')

#3.1

Iext_3 = function(t) {
  if (t < 0 || t > 700) return(0)
  
  if (t < 60) return(0)
  if (t < 300) return(0.6)
  if (t < 400) return(0)
  if ()
  
}


