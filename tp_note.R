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

plot(out[,2], type='l', ylim=c(-3, 3))
lines(out_euler[,2], type='l', col='red')
title("Van der Pol solution for x(t) ")
legend('topright', legend=c('ode', 'euler'),
       col=c('black', 'red'), lwd=2)


plot(out[,3], type='l', ylim=c(-3, 3))
lines(out_euler[,3], type='l', col='red')
title("Van der Pol solution for x'(t) ")
legend('topright', legend=c('ode', 'euler'),
       col=c('black', 'red'), lwd=2)

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
title("L'erreur entre euler et ode dans le temp T=50")
legend('topright', legend=c('ode', 'euler'),
       col=c('black', 'red'), lwd=2)


# 1.4 trouver un constante le lipschitz pour le f1 pour trouver le majoration d'erreur
yini = c(x=1, z=0)
times = seq(0, 50, 0.001)
system.time(out_euler <- euler(times=times, y=yini, func=eq_van_der_pol, parms = NULL))

# temp pour h = 10^(-3)  est 43s

F = function(x, z) {
  if (x >= 3 && x <= -3) {
    return(0)
  }
  if (z >= 3 && z <= -3) {
    return(0)
    }
  
  return (abs(-2*x*eps*w0*z - w0^2 + eps*w0*(1-x^2)))
}

a=seq(-3, 3, by=0.01)
b=seq(-3, 3, by =0.01)
r = outer(a, b, FUN=F)
L = max(r)

# constante de Lipschits 7.68
# on suppose que C est d ordre 1 
# on calcule le pas de la majoration de l'erreur globale 10^(-8)
hp = 10^(-8) / (1 *50*exp(7.68*50))
# pour avoir le temps dans seconds
tp = 10^(-3)/hp * 43
# ça donnera 10^(+175) des seconds (infiniment de temps)

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

out_ode <- ode(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)
out_2 <- euler(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)
plot(out_2[,2], type='l', col='red')
lines(out_ode[,2], col='black')
title("solution de Fitzhugn-Nagumo de v")
legend('topright', legend=c('ode', 'euler'),
       col=c('black', 'red'), lwd=2)

plot(out_2[,3], type='l', col='red')
lines(out_ode[,3], col='black')
title("solution de Fitzhugn-Nagumo de w")
legend('topright', legend=c('ode', 'euler'),
       col=c('black', 'red'), lwd=2)

#3.1

runge_kutta = function(times, y, func, parms) {
  state = y
  
  trace = matrix(c(times[1], state), nrow=1)
  
  for (i in seq(1, length(times)-1)) {
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
    
    trace = rbind(trace, c(times[i+1], c(state)))
  }
  trace
}

Iext_3 = function(t) {
  if (t < 0 || t > 700) return(0)
  
  if (t < 60) return(0)
  if (t < 300) return(0.6)
  if (t < 400) return(0)
  if (t < 600) return(0.5)
  if (t < 700) return(0) 
}

parms <- c(
  a = 0.95,
  b = 0.4,
  tau = 11,
  Iext = Iext_3
)
yini = c(v=0, w=0)
times = seq(0, 700, 1)

# discretization equidistante
out_3 <- runge_kutta(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)
out_3_euler <- euler(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)
out_3_ode <- ode(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)

plot(times, out_3[,2], type='l')
lines(times, out_3_euler[,2], type='l', col='red')
lines(times, out_3_ode[,2], type='l', col='blue')
title("solution de Fitzhugn-Nagumo de v avec condition de question 3 (h equidist)")
legend('topright', legend=c('runge_kutta', 'euler', 'ode'),
       col=c('black', 'red', 'blue'), lwd=2)

#q3.2 l erreur pour euler est infini, parce que c est la methode d orde 1 est la solution change tres rapide
# runge_kutta n a pas de probleme ici parce que c est la methode d ordre 4
err_ode_euler = max(out_3_euler - out_3_ode)
print(err_ode_euler)
err_ode_runge_kutta = max(out_3 - out_3_ode)
print(err_ode_runge_kutta)


# discretization specisalisé
times = c(seq(0, 60, 1), 
          seq(60, 350, 0.01),
          seq(350, 400, 1), 
          seq(400, 600, 0.1), 
          seq(600, 700, 1))

out_3 <- runge_kutta(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)
out_3_euler <- euler(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)
out_3_ode <- ode(times=times, y=yini, func=eq_fitzh_nagumo, parms = parms)

plot(times, out_3[,2], type='l')
lines(times, out_3_euler[,2], type='l', col='red')
lines(times, out_3_ode[,2], type='l', col='blue')
title("solution de Fitzhugn-Nagumo de v avec condition de question 3 (h speciale)")
legend('topright', legend=c('runge_kutta', 'euler', 'ode'),
       col=c('black', 'red', 'blue'), lwd=2)


#q3.2 avec discretization speciale on peut diminuer l'erreur est euler commence a marcher
# runge_kutta n a pas de probleme parce que c est la methode d ordre 4
err_ode_euler = max(out_3_euler - out_3_ode)
print(err_ode_euler)
err_ode_runge_kutta = max(out_3 - out_3_ode)
print(err_ode_runge_kutta)
