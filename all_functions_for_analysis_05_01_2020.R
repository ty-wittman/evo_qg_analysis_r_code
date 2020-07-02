
### functions written for analysis of genetic variance-covariance matricis including sex specific G matricies and those which include a between sex component, G with B. 

## computes the size of the genetic variance-covariance matrix
matsize = function(X){
  E = eigen(X)
  ev = E$values
  sev = sum(E$values)
  return(sev=sev)
}

### return the index of the sum of the eigen values of two matrices you are comparing
eigindex = function(X,Y){
  e1 = sum(eigen(X)$values)
  e2  = sum(eigen(Y)$values)
  ern = ifelse(e1>e2,((e1/e2)-1)*-1,(e2/e1)-1)
  return(ern)
}

## returns the difference in the sum of the eigen values of matrices you are comparing
eigdif = function(X,Y){
  e1 = sum(eigen(X)$values)
  e2  = sum(eigen(Y)$values)
  ed = e1-e2
  return(ed)
}

## computes random skewers analysis (Chevrud 1996) with random skewers drawn from a uniform distribution, provides more detailed output than RS_unif_quick. Put in two matrices you want to compare and the number of skewers. Can easily be vectorized to create a sampling distribuiton. 

RS_unif_full = function(X,Y,nsim){
  m = nrow(X)
  set.seed(493242)
  S <- matrix(runif(nsim * m, min = -2, max = 2), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  R1 <- matrix(rep(1,nsim*m), nsim, m)
  R2 <- matrix(rep(1,nsim*m), nsim, m)
  for (i in 1:nrow(R1)){
    R1[i,] = (S[i,])%*%X
  }
  for (i in 1:nrow(R2)){
    R2[i,] = (S[i,])%*%Y
  }
  vc = rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2))
  mvc = mean(rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2)))
  r1n = apply(R1,1,Norm)
  r2n = apply(R2,1,Norm)
  return(list(vc2 = vc2,vc=vc,mvc=mvc,r1n=r1n,r2n=r2n))
  
}


## quick version with only the mean vector correlation as the output. 

RS_unif_quick = function(X,Y,nsim){
  m = nrow(X)
  S <- matrix(runif(nsim * m, min = -2, max = 2), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  R1 <- matrix(rep(1,nsim*m), nsim, m)
  R2 <- matrix(rep(1,nsim*m), nsim, m)
  for (i in 1:nrow(R1)){
    R1[i,] = (S[i,])%*%X
  }
  for (i in 1:nrow(R2)){
    R2[i,] = (S[i,])%*%Y
  }
  mvc = mean(rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2)))
  return(mvc=mvc)
}

### compute the ratio of the magnitude of response to random skewers, can be used as a measure of the relative size of the matrices being compared. 

RS_unif_magnitude_ratio = function(X,Y,nsim){
  m = nrow(X)
  set.seed(4324352)
  S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  
  R1 = matrix(t(apply(S,1,"%*%",X)),ncol=m,nrow=nsim)
  R2 = matrix(t(apply(S,1,"%*%",Y)),ncol=m,nrow=nsim)
  
  r1n = apply(R1,1,Norm)
  r2n = apply(R2,1,Norm)
 meanrat = mean(r2n/r1n)
  return(meanrat)
}

### compute an index of the magnitude of response to random skewers, can be used as a measure of the relative size of the matrices being compared. 


RS_unif_mag_index = function(X,Y,nsim){
  m = nrow(X)
  set.seed(493242)
  S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  
  R1 = matrix(t(apply(S,1,"%*%",X)),ncol=m,nrow=nsim)
  R2 = matrix(t(apply(S,1,"%*%",Y)),ncol=m,nrow=nsim)
  
  r1n = apply(R1,1,Norm)
  r2n = apply(R2,1,Norm)
  rn = ifelse(r1n>r2n,(r1n/r2n)-1,((r2n/r1n)-1)*-1)
  mrn=mean(rn)
  return(mrn)
}


### Random skewers analysis but with skewers standardized to have the same selection intensity for both matrices, need to provide vectors of the phenotypic standard deviations for each trait for each matrix.


RS_std_intensity = function(X,Y,sdx,sdy,nsim){
  m = nrow(X)
  set.seed(493242)
  S <- matrix(runif(nsim * m, min = -2, max = 2), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  R1 <- matrix(rep(1,nsim*m), nsim, m)
  R2 <- matrix(rep(1,nsim*m), nsim, m)
  for (i in 1:nrow(R1)){
    R1[i,] = (S[i,]/sdx)%*%X
  }
  for (i in 1:nrow(R2)){
    R2[i,] = (S[i,]/sdy)%*%Y
  }
  mvc = mean(rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2)))
  return(mvc = mvc)
  
}

###  Same as above but also standardized response in terms of change in standard deviations of the mean

RS_std_intensity_and_response = function(X,Y,sdx,sdy,nsim){
  m = nrow(X)
  set.seed(493242)
  S <- matrix(runif(nsim * m, min = -2, max = 2), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  R1 <- matrix(rep(1,nsim*m), nsim, m)
  R2 <- matrix(rep(1,nsim*m), nsim, m)
  for (i in 1:nrow(R1)){
    R1[i,] = (S[i,]/sdx)%*%X
    R1[i,] = R1[i,]/sdx
  }
  for (i in 1:nrow(R2)){
    R2[i,] = (S[i,]/sdy)%*%Y
    R2[i,] = R2[i,]/sdy
  }
  mvc = mean(rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2)))
  return(mvc = mvc)
  
}

###  Random skewers are now treated as selectiton differentials and the phenotypic variance-covariance matrix is used to bend them into selection gradients. 



RS_std_intensity_gradients_from_P = function(X,Y,px,py,nsim){
  m = nrow(X)
  sdpx = sqrt(diag(px))
  sdpy = sqrt(diag(py))
  ipy = inv(py)
  ipx = inv(px)
  S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  Sx = matrix(t(apply(S,1,"*",sdpx)),ncol=m,nrow=nsim)
  Sxs = matrix(t(apply(Sx,1,"%*%",ipx)),ncol=m,nrow=nsim)
  
  Sy = matrix(t(apply(S,1,"*",sdpy)),ncol=m,nrow=nsim)
  Sys = matrix(t(apply(Sy,1,"%*%",ipy)),ncol=m,nrow=nsim)
  
  R1 = matrix(t(apply(Sxs,1,"%*%",X)),ncol=m,nrow=nsim)
  
  R2 = matrix(t(apply(Sys,1,"%*%",Y)),ncol=m,nrow=nsim)
  
  vc = rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2))
  
  mvc = mean(rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2)))
  return(mvc=mvc)
}

###  Same as above but response is standardized in terms of standard deviations


RS_std_intensity_and_response_gradients_from_P = function(X,Y,px,py,nsim){
  m = nrow(X)
  sdpx = sqrt(diag(px))
  sdpy = sqrt(diag(py))
  ipy = inv(py)
  ipx = inv(px)
  S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  Sx = matrix(t(apply(S,1,"*",sdpx)),ncol=m,nrow=nsim)
  Sxs = matrix(t(apply(Sx,1,"%*%",ipx)),ncol=m,nrow=nsim)
  
  Sy = matrix(t(apply(S,1,"*",sdpy)),ncol=m,nrow=nsim)
  Sys = matrix(t(apply(Sy,1,"%*%",ipy)),ncol=m,nrow=nsim)
  
  R1 = matrix(t(apply(Sxs,1,"%*%",X)),ncol=m,nrow=nsim)
  R11 = matrix(t(apply(R1,1,"/",sdpx)),ncol=m,nrow=nsim)
  R2 = matrix(t(apply(Sys,1,"%*%",Y)),ncol=m,nrow=nsim)
  R22 = matrix(t(apply(R2,1,"/",sdpy)),ncol=m,nrow=nsim)
  
  
  mvc = mean(rowMeans(R11*R22)/sqrt(rowMeans(R11^2)*rowMeans(R22^2)))
  return(mvc=mvc)
}

### function for the correlation between two vectors

VecCor = function(X,Y){
  return(rowMeans(X*Y)/sqrt(rowMeans(X^2)*rowMeans(Y^2)))
}

## random skewers analysis but for G with B, creates sexually antagonistic selection skewers and calculates response for males and females using both sex specific components G, and the between sex components B.
### So only a single matrix, which contains both male and female and between sex components, is imputed. 
## you can control the type of selection with the type argument, there are options for both sexually antagonistic and sexually concordant selection.
##OS, all random skewers for a shared trait are given the opposite sign, but are free to differ in magnitude
##OS-ALL, random skewers for a shared trait are opposite in sign and the same in magnitude, as sexually antagonisitic as you can get
##SS-DMM, same sign, double magnitude in males
##SS-DMF, same sign, double magnitude in females


SAskewer = function(X,nsim,type){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  set.seed(493242)
  S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  colnames(S) = colnames(X)
  R <- matrix(rep(1,nsim*m), nsim, m)
  colnames(R) = colnames(X)
  if(type == "OS-ALL"){
    for (i in seq(1,5)){
      S[,i]= ifelse(S[,1] >0 , abs(S[,i]),-abs(S[,i] ))
    }
    for (i in seq(1,5)){
      S[,i+5] = ifelse(S[,i] >0 , -abs(S[,i+5]),abs(S[,i+5] ))
    }
  } else {
    if(type=="OS"){
      for (i in seq(1,5)){
        S[,i+5]= ifelse(S[,i] >0 & S[,i+5] > 0 , -S[,i+5],S[,i+5] )
        S[,i+5] = ifelse(S[,i] <0 & S[,i+5] <0 , -S[,i+5],S[,i+5] )
      }
      
    } else {
      if(type=="SS-DMM"){
        for (i in seq(1,5)){
          S[i+5,]= ifelse(S[i,] >0  , S[i,]/2,S[i+5,] )
          S[i+5,] = ifelse(S[i,] <0 , S[i,]/2,S[i+5,] )
        }
        
        
      } else {
        if(type =="SS-DMF"){
          for (i in seq(1,5)){
            S[,i+5]= ifelse(S[,i] >0  , S[,i]*2,S[,i+5] )
            S[,i+5] = ifelse(S[,i] <0 , S[,i]*2,S[,i+5] )
          }
          
        } else {
          stop("enter type")
        }
      }
      
    }
  }
  for (i in 1:nrow(R)){
    R[i,] = 0.5*(S[i,])%*%X
  }
  Rtop = R[,1:(nrow(X)/2)]
  Rbot = R[,((nrow(X)/2)+1):nrow(X)]
  corbytrait = colMeans(Rtop*Rbot)/sqrt(colMeans(Rtop^2)*colMeans(Rbot^2))
  r = mean(corbytrait)
  vectorcordist = rowMeans(Rtop*Rbot)/sqrt(rowMeans(Rtop^2)*rowMeans(Rbot^2))
  meanvectorcor = mean(rowMeans(Rtop*Rbot)/(sqrt(rowMeans(Rtop^2))*sqrt(rowMeans(Rbot^2))))
  meanevol = mean(rowMeans(R*S)/(sqrt(rowMeans(R^2))*sqrt(rowMeans(S^2))))
  evodist = rowMeans(R*S)/(sqrt(rowMeans(R^2))*sqrt(rowMeans(S^2)))
  return(list(S = S,vcd = vectorcordist,mvc = meanvectorcor,r=r,corbytrait = corbytrait,Rtop = Rtop,Rbot = Rbot,evodist=evodist,meanevol=meanevol))
}


### faster version of the above code with limited output. 
SAskewer_quick = function(X,nsim,type){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  S <- matrix(runif(nsim * m, min = -2, max = 2), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  colnames(S) = colnames(X)
  R <- matrix(rep(1,nsim*m), nsim, m)
  colnames(R) = colnames(X)
  if(type == "OS-ALL"){
    for (i in seq(1,5)){
      S[,i]= ifelse(S[,1] >0 , abs(S[,i]),-abs(S[,i] ))
    }
    for (i in seq(1,5)){
      S[,i+5] = ifelse(S[,i] >0 , -abs(S[,i+5]),abs(S[,i+5] ))
    }
  } else {
    if(type=="OS"){
      for (i in seq(1,5)){
        S[,i+5]= ifelse(S[,i] >0 & S[,i+5] > 0 , -S[,i+5],S[,i+5] )
        S[,i+5] = ifelse(S[,i] <0 & S[,i+5] <0 , -S[,i+5],S[,i+5] )
      }
      
    } else {
      if(type=="SS-DMM"){
        for (i in seq(1,5)){
          S[i+5,]= ifelse(S[i,] >0  , S[i,]/2,S[i+5,] )
          S[i+5,] = ifelse(S[i,] <0 , S[i,]/2,S[i+5,] )
        }
        
        
      } else {
        if(type =="SS-DMF"){
          for (i in seq(1,5)){
            S[,i+5]= ifelse(S[,i] >0  , S[,i]*2,S[,i+5] )
            S[,i+5] = ifelse(S[,i] <0 , S[,i]*2,S[,i+5] )
          }
          
        } else {
          stop("enter type")
        }
      }
      
    }
  }
  for (i in 1:nrow(R)){
    R[i,] = 0.5*(S[i,])%*%X
  }
  Rtop = R[,1:(nrow(X)/2)]
  Rbot = R[,((nrow(X)/2)+1):nrow(X)]
  vectorcordist = rowMeans(Rtop*Rbot)/sqrt(rowMeans(Rtop^2)*rowMeans(Rbot^2))
  meanvectorcor = mean(rowMeans(Rtop*Rbot)/(sqrt(rowMeans(Rtop^2))*sqrt(rowMeans(Rbot^2))))

  return(meanvectorcor)
}


### performes the above analysis but for a single focal trait at a time in a background of random selection on the other traits in the matrix. Useful for investigating how the off-diagonal structure of the matrix can influence evolution of sexual dimorphism even when those traits are not under selection 
sadecomp = function(X,nsim,type){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  mat.names = paste("response_trait",1:m/2,sep="_")
  set.seed(56574946)
  Sd = array(runif(nsim * m, min = -1, max = 1),dim=c(nsim,m,m/2))
  tr = seq(1,m/2,1)
  trait_decomp = array(rep(1,nsim*m*m),dim=c(nsim,m,m/2))
  for(i in tr ){
    Sd[,i+5,i] = ifelse(Sd[,i,i] > 0 & Sd[,i+5,i] >0, -Sd[,i+5,i],Sd[,i+5,i])
    Sd[,i+5,i] = ifelse(Sd[,i,i] <0 & Sd[,i+5,i] <0 , -Sd[,i+5,i],Sd[,i+5,i] )
  }
  if(type == "uni"){
    for(i in tr){
      for(j in 1:10){
        for(x in 1:nsim){
          Sd[x,j,i] = ifelse(j==i|j-5 == i, Sd[x,j,i], 0)
        }
      }
    }
    for (j in tr){
      for (i in 1:nrow(trait_decomp[,,j])){
        trait_decomp[i,,j] = 0.5*Sd[i,,j]%*%X
      }
    }
  } else {
    for (j in tr){
      for (i in 1:nrow(trait_decomp[,,j])){
        trait_decomp[i,,j] = 0.5*Sd[i,,j]%*%X
      }
    }
  }
  trait_decomp <<- trait_decomp
  ssd <<- Sd
}



sa_decomp_uni = function(X,nsim){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  mat.names = paste("response_trait",1:m/2,sep="_")
  set.seed(564974956)
  Sd = array(runif(nsim * m, min = -1, max = 1),dim=c(nsim,m,m/2))
  tr = seq(1,m/2,1)
  trait_decomp = array(rep(1,nsim*m*m),dim=c(nsim,m,m/2))
  for(i in tr ){
    Sd[,i+5,i] = ifelse(Sd[,i,i] > 0 & Sd[,i+5,i] >0, -Sd[,i+5,i],Sd[,i+5,i])
    Sd[,i+5,i] = ifelse(Sd[,i,i] <0 & Sd[,i+5,i] <0 , -Sd[,i+5,i],Sd[,i+5,i] )
    
  }
  for (j in tr){
    for (i in 1:nrow(trait_decomp[,,j])){
      trait_decomp[i,,j] = 0.5*Sd[i,,j]%*%X
      
    }
    
  } 
  trait_decomp <<- trait_decomp
  Sd <<- Sd
}


