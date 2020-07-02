
### script load in all functions#####
matsize = function(X){
  E = eigen(X)
  ev = E$values
  sev = sum(E$values)
  return(sev=sev)
}

eigindex = function(X,Y){
  e1 = sum(eigen(X)$values)
  e2  = sum(eigen(Y)$values)
  ern = ifelse(e1>e2,((e1/e2)-1)*-1,(e2/e1)-1)
  return(ern)
}

eigdif = function(X,Y){
  e1 = sum(eigen(X)$values)
  e2  = sum(eigen(Y)$values)
  ed = e1-e2
  return(ed)
}


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



VecCor = function(X,Y){
  return(rowMeans(X*Y)/sqrt(rowMeans(X^2)*rowMeans(Y^2)))
}

selection_skewers = function(X,Y,px,py,nsim){
  m = nrow(X)
  sdpx = sqrt(diag(px))
  sdpy = sqrt(diag(py))
  ipy = inv(py)
  ipx = inv(px)
  S <- array(runif(nsim * m * m, min = -1, max = 1), dim=c(nsim, m,m))
  n = array(apply(S,3,function(S) sqrt(rowSums(S^2))),dim=c(nsim,1,m))
  for(i in 1:m){
    S[,i,i] = abs(S[,i,i])+0.2
  }
  for(i in 1:m){
    S[,,i]=S[,,i]/n[,,i]
  }
  Sx =  aperm(apply(S,c(1,3),"*",sdpx),c(2,1,3))
  Sxs = aperm(apply(Sx,c(1,3),"%*%",ipx),c(2,1,3))
  
  Sy =  aperm(apply(S,c(1,3),"*",sdpy),c(2,1,3))
  Sys = aperm(apply(Sy,c(1,3),"%*%",ipy),c(2,1,3))
  
  R1 = aperm(apply(Sxs,c(1,3),"%*%",X),c(2,1,3))
  R1 = aperm(apply(R1,c(1,3),"/",sdpx),c(2,1,3))
  
  R2 = aperm(apply(Sys,c(1,3),"%*%",Y),c(2,1,3))
  R2 = aperm(apply(R2,c(1,3),"/",sdpy),c(2,1,3))
  
  vc = sapply(1:m,function(i) VecCor(R1[,,i],R2[,,i]))
  
  return(list(vc=vc,Sxs=Sxs,Sys=Sys,R1=R1,R2=R2))
}


selection_skewers_uni = function(X,Y,px,py,nsim){
  m = nrow(X)
  sdpx = sqrt(diag(px))
  sdpy = sqrt(diag(py))
  ipy = inv(py)
  ipx = inv(px)
  S <- array(runif(nsim * m * m, min = -1, max = 1), dim=c(nsim, m,m))
  n = array(apply(S,3,function(S) sqrt(rowSums(S^2))),dim=c(nsim,1,m))
  for(i in 1:m){
    S[,i,i] = abs(S[,i,i])+0.2
    for (j in 1:m){
      if (j != i){
        S[,j,i] = S[,j,i]*0.5
      }
    }
  }
  for(i in 1:m){
    S[,,i]=S[,,i]/n[,,i]
  }
  Sx =  aperm(apply(S,c(1,3),"*",sdpx),c(2,1,3))
  Sxs = aperm(apply(Sx,c(1,3),"%*%",ipx),c(2,1,3))
  
  Sy =  aperm(apply(S,c(1,3),"*",sdpy),c(2,1,3))
  Sys = aperm(apply(Sy,c(1,3),"%*%",ipy),c(2,1,3))
  
  R1 = aperm(apply(Sxs,c(1,3),"%*%",X),c(2,1,3))
  R1 = aperm(apply(R1,c(1,3),"/",sdpx),c(2,1,3))
  
  R2 = aperm(apply(Sys,c(1,3),"%*%",Y),c(2,1,3))
  R2 = aperm(apply(R2,c(1,3),"/",sdpy),c(2,1,3))
  vc = sapply(1:m,function(i) VecCor(R1[,,i],R2[,,i]))
  
  return(list(vc=vc,Sxs=Sxs,Sys=Sys,R1=R1,R2=R2))
}




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


