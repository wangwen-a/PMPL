#### define SCAD####
SCAD = function(lambda,t){
  if(abs(t)<=lambda){
    a=lambda*abs(t)
  }
  if((abs(t)<=3.7*lambda)&(abs(t)>lambda)){
    a=-(t^2-2*3.7*lambda*abs(t)+lambda^2)/(2*(3.7-1))
  }
  if(abs(t)>3.7*lambda){
    a=0.5*(3.7+1)*lambda^2
  }
  return(a)
}



#### define BIC####
BIC=function(beta,y,x){
  n = nrow(x)
  residual0 = y-x%*%beta
  quant = quantile(residual0,probs=seq(0,1,0.25))
  h = 1.06*n^(-0.2)*min(sd(residual0),(quant[4]-quant[2])/1.34)
  #Silverman’s rule of thumb bandwidth for Gaussian kernel
  
  print(paste("h in BIC function:",h))
  
  a = numeric(n)
  for(i in 1:n){
    a[i] = sum( gaussionkernel((residual0-residual0[i])/h)/n/h )
  }
  
  df=sum(beta!=0) 
  BIC = -sum( log(a) )*2+log(n)*df*max(log(log(length(beta))),1)
  return(BIC)
}


#### define BIC with fixed h####
BIC_h=function(beta,y,x,h){
  n = nrow(x)
  residual0 = y-x%*%beta
  h=h
  
  
  a = numeric(n)
  for(i in 1:n){
    a[i] = sum( gaussionkernel((residual0-residual0[i])/h)/n/h )
  }
  
  df=sum(beta!=0) 
  BIC = -sum( log(a) )*2+log(n)*df*max(log(log(length(beta))),1)
  return(BIC)
}

#### define Gaussian kernel####
gaussionkernel=function(x)
{
  exp(-0.5*(x^2) )/ sqrt(2*pi)
}




#### define SCAD_MPL####
SCAD.MPL=function(y,X,beta0,lambda)
{
  #### estimation procedure
  n = nrow(X)
  residual0 = y-X%*%beta0
  quant = quantile(residual0,probs=seq(0,1,0.25))
  h = 1.06*n^(-0.2)*min(sd(residual0),(quant[4]-quant[2])/1.34)
  print(paste("h in SCAD.MPL function:",h))
  
  beta.selection1 = beta0
  beta.new = beta0
  beta.old = beta0
  
  a = numeric(n)
  penality = numeric(length(beta0))
  
  for (kkk in 1:15){
    for (k in 1:length(beta0)){
      beta.choices0 = seq(beta.selection1[k]-0.5,beta.selection1[k]+0.5,by=0.0625)
      results = numeric(0) 
      times =length(beta.choices0)+1  
      beta.choices = numeric(times) 
      beta.choices[1] = 0
      beta.choices[2:times] = beta.choices0  
      for(hh in 1:times){
        beta.new[k] = beta.choices[hh]
        residual.new = y-X%*%beta.new
        for(i in 1:n){
          a[i] = sum( gaussionkernel((residual.new-residual.new[i])/h)/n/h )
        }
        
        
        penality[k] = SCAD(lambda,beta.new[k])
        results[hh] = -sum( log(a) )+penality[k]*n
       
      }
      beta.new[k] = beta.choices[which.min(results)]

    }
    
    if (sum((beta.new-beta.old)^2)==0){
      break
    }
    beta.old = beta.new 
  }

  
  return(beta.new)
}


#### define SCAD_MPL with fixed h####
SCAD.MPL_h=function(y,X,beta0,lambda,h)
{

  n = nrow(X)
  residual0 = y-X%*%beta0
  h=h
  
  
  beta.selection1 = beta0
  beta.new = beta0
  beta.old = beta0
  
  a = numeric(n)
  penality = numeric(length(beta0))
  
  for (kkk in 1:15){
    for (k in 1:length(beta0)){
      beta.choices0 = seq(beta.selection1[k]-0.5,beta.selection1[k]+0.5,by=0.0625)
      results = numeric(0) 
      times =length(beta.choices0)+1  
      beta.choices = numeric(times) 
      beta.choices[1] = 0
      beta.choices[2:times] = beta.choices0  
      for(hh in 1:times){
        beta.new[k] = beta.choices[hh]
        residual.new = y-X%*%beta.new
        for(i in 1:n){
          a[i] = sum( gaussionkernel((residual.new-residual.new[i])/h)/n/h )
        }
  
        penality[k] = SCAD(lambda,beta.new[k])
        results[hh] = -sum( log(a) )+penality[k]*n
        
      }
      beta.new[k] = beta.choices[which.min(results)]
      
    }
    
    if (sum((beta.new-beta.old)^2)==0){
      break
    }
    beta.old = beta.new 
  }
  
  
  return(beta.new)
}
        
        
#### define SCAD_MPL_BIC####
SCAD.MPL.BIC = function(y,X,beta0,lam.start,lam.end,lam.range){
  lambda.choices = seq(lam.start,lam.end,by=lam.range) 
  num.lambda = length(lambda.choices)
  BIC.results = numeric(0)
  for (num in 1:num.lambda){
    lambda = lambda.choices[num]
    beta.new = SCAD.MPL(y,X,beta0,lambda)
    
    
    BIC.results[num] = BIC(beta.new,y,X)
  }
  lambda.final = lambda.choices[which.min(BIC.results)]
  
  
  beta.final = SCAD.MPL(y,X,beta0,lambda.final)
  list=numeric(0)
  list[1]=lambda.final
  list[2:(length(beta0)+1)]=beta.final
  
  return(list)
}
### out:
## list[1]: the value of lambda
## list[2]: the value of beta estimated

###### define SCAD_MPL_BIC with fixed h####
SCAD.MPL.BIC_h = function(y,X,beta0,lam.start,lam.end,lam.range,h){
  lambda.choices = seq(lam.start,lam.end,by=lam.range)
  num.lambda = length(lambda.choices)
  BIC.results = numeric(0)
  for (num in 1:num.lambda){
    lambda = lambda.choices[num]
    beta.new = SCAD.MPL_h(y,X,beta0,lambda,h)
    
    BIC.results[num] = BIC_h(beta.new,y,X,h)
  }
  lambda.final = lambda.choices[which.min(BIC.results)]
  
  beta.final = SCAD.MPL_h(y,X,beta0,lambda.final,h)
  list=numeric(0)
  list[1]=lambda.final
  list[2:(length(beta0)+1)]=beta.final
  
  return(list)
}
### out:
## list[1]: the value of lambda
## list[2]: the value of beta estimated


