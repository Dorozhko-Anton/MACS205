##################################
######## Exercice 1 ##############
######## Méthode des trapezes ####

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

find_order(trapezeInt)


##################################
######## Exercice 2 ##############
######## Division du pas  des trapezes par deux ####


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

### TEST
p4 = function(x){x^4}
M =4
myfun = p4 
Qh = trapezeInt(myfun, 0, 1, M)
refineQh = refineTrapeze(myfun,0,1,M, Qh)
Qh2 = trapezeInt(myfun, 0, 1, 2*M)
err = Qh2 - refineQh
err


##################################
######## Exercice 4 ##############
######## Des trapèzes à Simpson ##
##################################


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

## TODO: coder le code pour tester la fonction dont on connait l'integral
test = function(FUN, a, b, realInt) {

  MM = 5:20
  Resultats = rep(0,length(MM));
  for (i in  1:length(MM)){
    Resultats[i]  = (b-a)*simpsonInt(FUN, a, b,MM[i]);
  }
  plot(MM,trueInt-Resultats);
  title("error as a function of M")
} 

## TEST
d = 4
a = 1; b = 5;

MyPolynome = function(x){x^d}
trueInt = (b^(d+1) - a^(d+1))/(d+1);

test(MyPolynome, a, b, trueInt)



## TEST
d = 4
a = 1; b = 5;
# a = pi/2; b=2*pi + pi/2

MyPolynome = function(x){x^d}
trueInt = (b^(d+1) - a^(d+1))/(d+1);

MM = 5:20
Resultats = rep(0,length(MM));
for (i in  1:length(MM)){
  Resultats[i]  = simpsonInt(MyPolynome,a,b,MM[i]);
  # Resultats[i]  = simpsonInt(MyPolynome,pi/2,2*pi+pi/2,MM(i));
}

plot(MM,trueInt-Resultats);
title("error as a function of M")

if(max(abs(trueInt-Resultats)) > trueInt * .Machine$double.eps*10 ){
  ## la condition satisfaite à partir d=4
  logerr = log(abs(trueInt-Resultats));
  plot(log(MM),logerr);
  title("log error as  a function of log(M)")
  neval = length(logerr)
  pente = (logerr[neval] - logerr[1])/(log(MM[neval]) - log(MM[1]))
  pente
}


##################################################
############# Exercice 5: evaluation de l'erreur a posteriori
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

## test
d = 13.4 ; M=9
a = 0.5; b=3;
pol = function(x){x^d}
trueInt = (b^(d+1) - a^(d+1))/(d+1)
vv= evalErrSimpson(pol,a,b, M)
q = vv[2] ; E = vv[1]
estErr = E
trueErr = q - trueInt 
relativeErrorOnError = (trueErr - estErr)/trueErr
relativeErrorOnError


########## Exercice 6: Richardson
############################

richardson = function(FUN,n,t,delta){
  ## Calcule le tableau des differences  divisees en 0 du 
  ## polynome d'interpolation en t,delta t, ... delta^n t
  ## renvoie un vecteur de taille n+1:
  ## le vecteur des A_{k,k}, k= 0 .. n 
  ## (pas la matrice).   
  ## La meilleure approximation est le dernier element A[n+1].
  ##
  lx = log(t)  +  log(delta) *(0:n)
  x = exp(lx) 
  A = sapply(x,FUN) 
  for( j in 2:(n+1)){
    A[j : (n+1) ] =  (A[j:(n+1)] - delta^(j-1)* A[(j-1):n])/ (1 - delta^(j-1))
  }
  return(A)
}

## test
myfun = function(x){sin(x)+ (cos(x)-1)*sin(2*x)}
n =10 ; t =pi/4 ;  delta = 1/4
A = richardson(myfun,n,t,delta)
lx = log(t) +  log(delta)*(0:n)
x = exp(lx);
y = sapply(x,myfun);

plot(0:n, y,col='blue', type = "l")
lines( 0:n,A,col='red');
legend('topright', legend=c('erreur naive', 'erreur Richardson'),
       col=c('blue', 'red'), lwd=2)

dev.new()
LerrRich = log(abs(A - myfun(0))) ;  
LerrNaive = log(abs(y - myfun(0)));
plot(-lx, LerrNaive,col='blue',type='l')
lines(-lx, LerrRich,col='red')
grid()
legend('topright', legend=c('log-erreur naive','log-erreur Richardson'),
       col=c('blue', 'red'), lwd=2)



##########################
######## exercice 7:  Romberg


romberg =function(FUN,n,a,b,M){## methode de Romberg avec n etapes
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


## test que Simposon de pas h est equivalent à trapèze de pas h + Romberg d'ordre de rang 1
myfun = function(x) {cos(x)}
n = 1; M = 3; a = 0; b = pi/2
I1 = romberg(myfun, n, a, b, M)
I2 = simpsonInt(myfun, a, b, M)
diff = I1[n+1] - I2
diff

## test
d=6
myfun = function(x){x^d}
n =10; M = 5; a= 0 ; b = 1.5
A = romberg(myfun, n, a, b , M)
B = rep(0,n+1)
B[1] = trapezeInt(myfun,a,b,M) 
Mc = M
for( i in  2:(n+1)){
  B[i] = refineTrapeze( myfun,a,b, Mc, B[i-1])
  Mc = 2* Mc 
}

plot(0:n, B,col='blue', type='l')
lines( 0:n, A,col = 'red');

LerrRich = log(abs(A - (b^(d+1) - a^(d+1))/(d+1))) ;  
LerrNaive = log(abs(B - (b^(d+1) - a^(d+1))/(d+1)));
plot((0:n), LerrNaive/log(4),col='blue',type='l')
lines( (0:n), LerrRich/log(4), col='red')
grid()
legend('topright', legend=c('log-erreur naive', 'log-erreur Romberg'),
       col=c('blue', 'red'), lwd=2)
