
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
r = dividif(c(1, 2, 3), c(1, 2, 3))
stopifnot(r == c(1, 1, 0))


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
  f  = a[n] * (z-x[n-1]) + a[n-1]
  for( i in 2:(n-1)){
    ## Complete the loop
    f = f*(z-x[n-i]) + a[n-i]
  }
  return(f)
}
 
# 
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

n = 10
x = seq(-1, 1, length.out = n)
myfun = function(x) {
  cos(5*x/(2 + x))* 1/(1 + (1*x)^2) 
}
y = sapply(x, myfun)

x_dense = seq(-1, 1, length.out = 10*n)
y_dense = sapply(x_dense, myfun)

plot(x_dense, y_dense, type="l")
lines(x_dense, interpolDividif(x, y, x_dense), col='red')

#ex5
funRunge = function(x) {
  1/(1+(10*x)^2)
}

n = 40
x = seq(-1, 1, length.out = n)
y = sapply(x, funRunge)

x_dense = seq(-1, 1, length.out = 10*n)
y_dense = sapply(x_dense, funRunge)

plot(x_dense, y_dense, type="l")
lines(x_dense, interpolDividif(x, y, x_dense), col='red')


# Tchebyshev nodes
genTchebyshevNodes = function(n) {
  cos(((0:(n-1) + 0.5)*pi/n))
}
n = 20
x = genTchebyshevNodes(n+1)
y = sapply(x, funRunge)

x_dense = seq(-1, 1, length.out = 10*n)
y_dense = sapply(x_dense, funRunge)

plot(x_dense, y_dense, type="l", col='blue')
lines(x_dense, interpolDividif(x, y, x_dense), col='red')


# 
# 
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
interpolLagrange(n=10, a=-1, b=2, neval=80,  nodes = 'cheby', FUN=funRunge, Plot=TRUE)

# 
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
piecewiseInterpol(n=10, nInt=5, a=-1, b=2, neval=80,  nodes = 'cheby', FUN=funRunge, Plot=TRUE)


interpolLagrange(n=60, a=-1, b=2, neval=80,  nodes = 'cheby', FUN=funRunge, Plot=TRUE)
piecewiseInterpol(n=4, nInt=15, a=-1, b=2, neval=80,  nodes = 'cheby', FUN=funRunge, Plot=TRUE)
