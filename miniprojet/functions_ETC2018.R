space_int <- 2^(-8)
load('hist1.rda')
apar <- 1.7 ## parametre alpha de la loi beta
bpar <- 5.1 ## parametre alpha de la loi beta

hist_value <- function(histog, x)
{
  if(any(x<0)|| any( x>=1))
  {
    stop('x n est pas dans [0,1]')
  }
  else{
    index_x <- floor(x/space_int);
    y <- histog[index_x+1];
    return(y)
  }
}


dividif=function(x,y){
  ##  Newton's Divided differences
  ## @param x: a vector containing the interpolation nodes 
  ## @param y: a vector of same size as x:
  ##           values of the interpolated function at the nodes
  ## @return : a vector of same size as x:
  ##          the divided differences
  ##          \eqn{f_[x_0, ... x_k]} of order 'length(x) -1'. 
  
  n = length(x) - 1
  d  = y 
  for (j in 2:(n+1) ) {
    d[j : (n+1) ] = ( d[j:(n+1)] - d[(j-1):n])/(x[j:(n+1)] - x[1:(n-j+2)])
  }
  return(d)
}

hornerNewton = function(a,x,z){
  ## Horner's method: Evaluates  a polynom P at points z, given
  ## nodes x and the coefficients a of P in Newton's basis
  ##
  ## @param a : vector: the  coefficients of the polynomial in
  ##           Newton's basis
  ## @param x : the interpolation nodes.
  ## @param z : vector of points where the polynom needs to be
  ##            evaluated.
  ## @return  : a vector of same size as z: the value of the
  ##            polynomial at points z.
  ##
  n = length(x);
  f  = a[n] * rep(1, length(z))
  if (n==0) {stop('at least one interpolating point is needed')}
  if (n >= 2) {
    for( i in 1:(n-1)){
      ## Complete the loop
      f = f*(z-x[n-i]) + a[n-i]
    }
  }
  return(f)
}

interpolDividif=function(x,y,z){
  ## Efficient Lagrange interpolation using Horner's method with
  ## Newton basis for evaluation
  ## @param x : vector containing the interpolation nodes
  ## @param y : vector of same size as x: values of the interpolated
  ##            function at the nodes
  ## @param z : vector of points where the  interpolating polynomial
  ##            needs to be evaluated.
  ## @return  : vector of same size as z: the value of the
  ##            interpolating polynomial at points z.
  
  ## Complete the code
  coefs = dividif(x, y)
  return (hornerNewton(coefs, x, z))
}

# Tchebyshev nodes
genTchebyshevNodes = function(n) {
  cos(((0:(n-1) + 0.5)*pi/n))
}

interpolLagrange =function(n, a, b, neval, nodes = 'equi', FUN, Plot){
  ## Generic Lagrange interpolation, with equidistant or Chebyshev nodes.
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param a : left end-point of the interval
  ## @param b : right end-point of the interval
  ## @param neval :number of evaluation points (a regular grid will be
  ## used on [a,b]
  ## @param nodes :string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## using Chebyshev nodes.
  ## @param FUN: the function to be interpolated
  ## @param Plot : logical. Setting 'Plot' to TRUE produces a plot
  ## showing the graph of
  ## the true functions and its interpolation.
  ## @return : vector of size 'neval': the values of the Lagrange
  ## polynomial on an equi-distant grid.
  
  if (nodes == "equi"){
    x =  seq(-1, 1, length.out = n+1)
  }
  else if (nodes == "cheby"){
    x = genTchebyshevNodes(n+1)
  }
  else{stop("the nodes must be either 'equi' or 'cheby'") }
  
  ##
  ## Complete the code: compute a vector 'f' containing
  ## the interpolated  values on an equidistant
  ## evaluation grid 'z'.
  ##
  ##
  z = seq(a, b, length.out = neval)
  u = (z - (a+b)/2)*2/(b-a)
  y = sapply((a+b)/2 + (b-a)/2*x, FUN = FUN)
  f = interpolDividif(x, y, u)
  
  if( Plot ){
    if (nodes == "equi"){ methodName = " equidistant "}
    else {   methodName = " Chebyshev "}
    
    plot(z, sapply(z,FUN), type="l", ylim=range(c(y,f)), col='black')
    title(main = paste("Lagrange interpolation with ",
                       toString(n+1), methodName,
                       " nodes", sep=""))
    lines(z,f, col = 'red')
    
    legend('topright', legend=c('function','interpolation'),
           col = c('black','red'), lwd=1)
    
  }
  return(f)
}

piecewiseInterpol=function(n,nInt,a,b,neval, nodes = "equi", FUN, Plot){
  ## @param n : the degree of the interpolating polynomial on each
  ## subinterval
  ## @param nInt :  the number of sub-intervals
  ## @param a, b : endpoints of the interval
  ## @param neval : the number of points on the interpolating grid (on
  ## each subinterval)
  ## @param nodes : string, either "equi" (default) for equidistant
  ## Lagrange interpolation (on each subinterval) or "cheby" for
  ## chebyshev nodes.
  ## @param FUN the function to be interpolated
  ## @param Plot : logical. Should the result be plotted ?
  ## @return : a matrix with 2 rows and neval * nInt -neval + 1:
  ## values of the interpolated funtion on a regular grid (first row)
  ## and the corresponding abscissas (second row).
  
  intEndPoints = seq(a,b,length.out = nInt+1)
  f = c()
  z = c()
  for (m in 1:nInt){
    A = intEndPoints[m]; B = intEndPoints[m+1]
    
    fm = interpolLagrange(n=n, a=A, b=B, neval=neval,  nodes = nodes, FUN=FUN, Plot=FALSE)
    zm = seq(A, B, length.out = neval)
    
    if( m >= 2){
      ## remove first element of zm, fm to avoid
      ## duplicate values of the  interpolating vector
      fm = fm[-1]
      zm = zm[-1]
    }
    z = c(z,zm)
    f = c(f,fm)
  }
  
  if (Plot){
    if (nodes == "equi") {methodName = " equidistant "}
    else  {methodName = " Chebyshev "}
    
    
    plot(z, sapply(z,FUN),type="l")
    title(main = paste("Piecewise  Lagrange  interpolation with ",
                       toString(n+1), methodName, " nodes  on ",
                       toString(nInt), " Intervals", sep=""))
    lines(z,f, col='red', lwd=2)
    legend('topright', legend = c('function','interpolation'),
           lwd=c(1,2), col=c('black','red'))
  }
  return(rbind(f,z) )
}


# 2.1
densite = function (x) {
  phi = dbeta(x, shape1 = 1.7, shape2 = 5.1)
  return (phi)
}

# Test densite
# x = seq(0,1,0.001)
# f = densite(x)
# plot(x, f)

f_hist_value = function(x) {
  hist_value(histog = hist1, x=x)
}


## Méthodes de quadrature

trapezeInt =function(FUN,a,b,M){
  ##' TRAPEZOIDAL INTEGRATION RULE (COMPOSITE)
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : number of intervals (each of size (b-a)/M)
  ##' @return: the value of the composite trapezoidal quadrature. 
  x = seq(a,b, length.out= M+1)
  y = sapply(x, FUN)
  h = (b-a)/M 
  q = h*(-0.5*y[1] + sum(y) - 0.5*y[length(y)])
  
  return(q)
}

find_order = function(FUN, tolerance=1e-7) {
  p = 0
  while (TRUE) {
    f = function (x) x**p
    int_lib = integrate(f, lower=-1, upper=1)
    int_to_check = FUN(f, -1, 1, 1)
    
    if (abs(int_to_check - int_lib$value) > tolerance) {
      break
    }
    p = p + 1
  }
  return (p-1)
}

refineTrapeze=function(FUN,a,b,M,q){
  ##' refinement of the subdivision step: incremental method
  ##' @param FUN : the function to be integrated
  ##' @param a, b : interval end points 
  ##' @param M : initial number of intervals (each of size (b-a)/M)
  ##'  having been used to compute q
  ##' @param  q : the value of the trapezoidal  quadrature method
  ##'  of stepsize (b-a)/M
  ##' @return : the value of the quadrature for a stepsize h' = h/2
  h = (b-a)/M
  ##  x : a vector of size M :
  ##     the additional abscissas where 'fun' must be evaluated.
  x = a + seq(1, 2*M-1, by=2)*h/2   
  y = sapply(x, FUN)
  Q = 0.5*(q + h*sum(y))
  return(Q)
}

simpsonInt = function(FUN,a,b,M){
  ##' Simpson integration via trapeze rule
  ##' uses the fact that 
  ##' simpson(h) = 4/3(trapeze(h/2) - 1/4 trapeze(h))
  h = (b-a)/M;
  qtrapeze = trapezeInt(FUN, a, b, M)  
  qrefined = refineTrapeze(FUN, a, b, M, qtrapeze) 
  q =  4/3 * (qrefined - 1/4 * qtrapeze) 
  return(q)
}

evalErrSimpson=function(FUN,a,b,M){
  ## Computes an approximation E of the error 
  ## for the composite Simpson rule of step h=(b-a)/(2M). 
  ##This requires computing I_M and I_{2M}. 
  ##The value  q = I_{2M} is also returned. 
  qth = trapezeInt(FUN,a,b,M)   ## M +1 evaluations
  qth2 = refineTrapeze ( FUN,a,b,M,qth )  ## M evaluations
  qth4 = refineTrapeze ( FUN,a,b,2*M,qth2 )   ## 2M evaluations
  simps_h =   4/3*(qth2 - 1/4* qth ) 
  simps_h2 =  4/3*(qth4 - 1/4* qth2 ) 
  q = simps_h2  
  E = (simps_h - simps_h2)/15 
  return(c(E,q))
}


chooseM = function(FUN, a, b, tolerance) {
  ##' estimate M - number of intervals to interpolate with given tolerance
  ##' @param FUN : funciton to interpolate
  ##' @param a : start of interval of interpolation
  ##' @param b : end of interval of interpolation
  ##' @param tolerance : absolute error
  ##' @return M : return number of intervals
  M = 10
  error = 1
  Mvector = c()
  Qvector = c()
  Evector = c()
  
  while (2*error > tolerance) {
    res = evalErrSimpson(FUN, a, b, M)
    
    Mvector = c(Mvector, M)
    Qvector = c(Qvector, res[2])
    error = abs(res[1])
    Evector = c(Evector, error)
    M = 2*M
  }
  
  return (rbind(Mvector, Qvector, Evector))
}

densite_bruitee = function(x) {
  hist_value(histog = hist1, x=x)
}

romberg =function(FUN,n,a,b,M){
  ## methode de Romberg avec n etapes
  ## appliquee sur la fonction FUN sur l'intervalle (a,b), avec un
  ## pas initial h = (b-a)/M
  h= (b-a)/M 
  A = rep(0, n+1)
  A[1] = trapezeInt(FUN,a,b,M);
  Mc = M
  ## initialisation des differences divisees
  for( i in 2:(n+1)){
    A[i] = refineTrapeze( FUN,a,b, Mc, q= A[i-1])
    Mc = 2*Mc 
  }
  delta = 1/4;
  for (j in 2:(n+1)){
    A[j : (n+1) ] = (A[j:(n+1)] - delta^(j-1)* A[(j-1):n])/ (1 - delta^(j-1))
  }
  return(A)
}

