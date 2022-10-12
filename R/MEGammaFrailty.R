## generate sample # TWO sample?

nsam <- function(a1,a2,th,be,n)
{
  u = rgamma(n, 1/th, rate=1/th)
  
  q = length(be)
  
  x = array( runif(2*n*q, min=0, max=0.5), dim=c(2,n,q) )
  
  Be = matrix(rep(be,each=n),n,q)
  
  T = matrix(0,2,n)
  
  cen = runif(n,0,3)    ## (0,3)/100=10%  
  
  
  U1 = runif(n,0,1)    
  U2 = runif(n,0,1)
  
  T[1,] <- -log(U1)/(a1*u*exp(apply(x[1,,]*Be,1,sum))) #?
  T[2,] <- ( exp( -log(U2)/(u*exp(apply(x[2,,]*Be,1,sum))) )-1 )/a2 #?
  
  d = 1*(T<=cen)
  y=pmin(T,cen)
  
  la1 = a1*y[1,]
  la2 = log(1+a2*y[2,])
  
  return(list(y=y,d=d,x=x,la1=la1,la2=la2,w=u))
}


############### log-likelihood function ############


loggamma = function(x,d,y,la1,la2,th,be,n)
{
  La1 = La2 = rep(0,n)
  for(i in 1:n){
    La1[i] = sum( la1*(y[1,]<=y[1,i]) )
    La2[i] = sum( la2*(y[2,]<=y[2,i]) )
  }
  
  q = length(be)
  
  BE = matrix(rep(be,each=n),n,q)
  C = 1/th + La1*exp(rowSums(x[1,,]*BE)) + La2*exp(rowSums(x[2,,]*BE))
  
  A = 1/th+colSums(d)
  AC = A/C
  
  d1 = d[1,]
  d2 = d[2,]
  
  l1 = sum(lgamma(A))-n*(lgamma(1/th)+log(th)/th)- sum(A*log(C)) 
  l2 = sum( log(la1[d1!=0]) ) + sum( log(la2[d2!=0]) ) 
  l3 = sum( d[1,]*rowSums(x[1,,]*BE) + d[2,]*rowSums(x[2,,]*BE) )
  
  ell = l1+l2+l3
  
  
  return(ell)
  
}


###################  F function  ###################################


F = function(la1,la2,be,th){
  
  La1 = La2 = rep(0,n)
  for(i in 1:n){
    La1[i] = sum( la1*(y[1,]<=y[1,i]) )
    La2[i] = sum( la2*(y[2,]<=y[2,i]) )
  }
  
  # la1-la2 #
  
  BE = matrix(rep(be,each=n),n,q)
  
  C0 = La1*exp(rowSums(x[1,,]*BE)) + La2*exp(rowSums(x[2,,]*BE))
  
  C = 1/th + La1*exp(rowSums(x[1,,]*BE)) + La2*exp(rowSums(x[2,,]*BE))
  
  A = 1/th+colSums(d)
  AC = A/C
  
  
  E_01 = AC*exp(rowSums(x[1,,]*BE)) 
  SUM_01 = cumsum( (E_01[order(y[1,])])[seq(n,1,-1)] )
  SUM_01 = ( SUM_01[seq(n,1,-1)] )[rank(y[1,])] 
  
  la1 = d[1,]/SUM_01
  if(any(is.na(la1))) la1 = la1.ini
  
  
  E_02 = AC*exp(rowSums(x[2,,]*BE)) 
  SUM_02 = cumsum( (E_02[order(y[2,])])[seq(n,1,-1)] )
  SUM_02 = ( SUM_02[seq(n,1,-1)] )[rank(y[2,])] 
  
  la2 = d[2,]/SUM_02
  if(any(is.na(la2))) la2 = la2.ini
  
  
  
  be.est = rep(0,q)
  
  for(p in 1:q){
    
    E1 =  La1*x[1,,p]*exp(rowSums(x[1,,]*BE)) +  La2*x[2,,p]*exp(rowSums(x[2,,]*BE))
    DE_1 = sum( d*x[,,p] ) - sum( AC*E1 )  
    
    AVE_X = apply(abs(x),c(1,2),sum)/abs(x[,,p])
    E2 = La1*x[1,,p]^2*AVE_X[1,]*exp(rowSums(x[1,,]*BE)) + La2*x[2,,p]^2*AVE_X[2,]*exp(rowSums(x[2,,]*BE))
    DE_2 = -2*sum( (colSums(d)+2/th)*E2/C )   
    
    be.est[p] = be[p] - DE_1/DE_2
    
  }    
  
  be = be.est 
  
  # th #
  Q01 = n*(digamma(1/th)+log(th)-2)/(th^2) + sum( log(C)+colSums(d)/C-digamma(A)+2/(C*th) )/(th^2) + sum( C0/C )/(th^2)
  
  Q02 = n*(4-2*digamma(1/th)-trigamma(1/th)/th-log(th))/(th^3)+2*sum(trigamma(A)/(2*th)+digamma(A)-log(C)-colSums(d)/C-3/(C*th))/(th^3) - 5*sum( C0/C )/(th^3)
  
  th = th - Q01/Q02
  if(th>0) { th = th }
  
  return(list(la1=la1,la2=la2,be=be,th=th))
  
}


##############

N = 500

mm.gamma = matrix(0,N,23)

for(j in 1:N){
  
  be0 = c(rep(-2,10), rep(3,10))
  q = length(be0)
  n = 100
  
  yy = nsam(3,5,1,be0,n)
  y = yy$y 
  d = yy$d
  x = yy$x
  
  dim(x)
  sum(d[1,])/n
  sum(d[2,])/n
  
  
  ###################  initial value 
  
  
  la1 = rep(1/n,n)
  la2 = rep(1/n,n)
  th = 0.5
  be = rep(0.5,q)
  
  ###  
  
  ell=rep(0,10000)
  k=1
  
  ell[k] = loggamma(x,d,y,la1,la2,n,be,n)
  
  diff = 3
  start = proc.time()[1]
  
  while( diff > 0.0000001 ) 
  {
    Fa = F(la1,la2,be,th)
    la1 = Fa$la1
    la2 = Fa$la2
    be1 = Fa$be
    th = Fa$th
    
    u_be = be1 - be
    
    FFa = F(la1,la2,be1,th)
    la1 = FFa$la1
    la2 = FFa$la2
    be2 = FFa$be
    th = FFa$th
    
    v_be = be2 - 2*be1 + be
    
    al_be = sum(u_be*v_be)/sum(v_be^2)
    if(al_be>-1)
    {al_be = -1}
    
    be = be2 -2*al_be*u_be + al_be^2*v_be
    
    
    La1 = La2 = rep(0,n)
    for(i in 1:n){
      La1[i] = sum( la1*(y[1,]<=y[1,i]) )
      La2[i] = sum( la2*(y[2,]<=y[2,i]) )
    }
    
    # la1-la2 #
    
    BE = matrix(rep(be,each=n),n,q)
    
    C0 = La1*exp(rowSums(x[1,,]*BE)) + La2*exp(rowSums(x[2,,]*BE))
    
    C = 1/th + La1*exp(rowSums(x[1,,]*BE)) + La2*exp(rowSums(x[2,,]*BE))
    
    A = 1/th+colSums(d)
    AC = A/C
    
    
    E_01 = AC*exp(rowSums(x[1,,]*BE)) 
    SUM_01 = cumsum( (E_01[order(y[1,])])[seq(n,1,-1)] )
    SUM_01 = ( SUM_01[seq(n,1,-1)] )[rank(y[1,])] 
    
    la1 = d[1,]/SUM_01
    if(any(is.na(la1))) la1 = la1.ini
    
    
    E_02 = AC*exp(rowSums(x[2,,]*BE)) 
    SUM_02 = cumsum( (E_02[order(y[2,])])[seq(n,1,-1)] )
    SUM_02 = ( SUM_02[seq(n,1,-1)] )[rank(y[2,])] 
    
    la2 = d[2,]/SUM_02
    if(any(is.na(la2))) la2 = la2.ini
    
    
    
    # th #
    Q01 = n*(digamma(1/th)+log(th)-2)/(th^2) + sum( log(C)+colSums(d)/C-digamma(A)+2/(C*th) )/(th^2) + sum( C0/C )/(th^2)
    
    Q02 = n*(4-2*digamma(1/th)-trigamma(1/th)/th-log(th))/(th^3)+2*sum(trigamma(A)/(2*th)+digamma(A)-log(C)-colSums(d)/C-3/(C*th))/(th^3) - 5*sum( C0/C )/(th^3)
    
    th1 = th - Q01/Q02
    if(th1>0) { th = th1 }
    
    
    
    ell[k+1]=loggamma(x,d,y,la1,la2,n,be,n)
    
    diff = abs(ell[k+1]-ell[k])/(1+abs(ell[k]))
    
    k = k+1
  }
  
  end = proc.time()[1]
  time = end-start 
  
  c(be,th,time)
  
  ell.est = loggamma(x,d,y,la1,la2,th,be,n)
  
  mm.gamma[j,] = c(be,th,time,ell.est)     # length = 23
  
}



Time = apply(mm.gamma,2,mean)[22]
Time

## MLE ##
apply(mm.gamma,2,mean)[c(1,5,10,15,20,21)]


## BIAS ##
abs( apply(mm.gamma,2,mean)[c(1,5,10,15,20,21)] - c(rep(-2,3), rep(3,2),1) )


## SD ##
apply(mm.gamma,2,sd)[c(1,5,10,15,20,21)]


