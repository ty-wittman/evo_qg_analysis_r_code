
### script load in all functions#####
### ave magnitude of between sex genetic correlation

rmf = function(X,Y){
  S = mean(diag(X[6:10,1:5])-diag(Y[6:10,1:5]))
  return(S)
}

rmf_u = function(X,Y){
  S = mean(diag(X[6:10,1:5]))
  return(S)
}

T_comp = function(X,Y){
  xod = upperTriangle(X,diag=F,byrow=T)
  yod = upperTriangle(Y,diag=F,byrow=T)
  S = abs(xod-yod)
  return(S=S)
}

#### t comp magnitude of correlations
### mean mag
mm = function(X){
  xod = upperTriangle(X,diag=F,byrow=T)
  m = sum(abs(xod))
  return(m)
}


T_comp_mag = function(X,Y){
  xod = upperTriangle(X,diag=F,byrow=T)
  yod = upperTriangle(Y,diag=F,byrow=T)
  S = abs(xod)-abs(yod)
  return(S=S)
}


T_comp_b = function(X,Y){
  xod = upperTriangle(X,diag=F,byrow=T)
  xld = lowerTriangle(X,diag=F,byrow=F)
  x = c(xod,xld)
  yod = upperTriangle(Y,diag=F,byrow=T)
  yld = lowerTriangle(Y,diag=F,byrow=F)
  y = c(yod,yld)
  S = mean(abs(x-y))
  return(S=S)
}

T_comp_b_diag = function(X,Y){
  xod = upperTriangle(X,diag=T,byrow=T)
  xld = lowerTriangle(X,diag=F,byrow=F)
  x = c(xod,xld)
  yod = upperTriangle(Y,diag=T,byrow=T)
  yld = lowerTriangle(Y,diag=F,byrow=F)
  y = c(yod,yld)
  S = mean(abs(x-y))
  return(S=S)
}

T_comp_b2 = function(X){
  xod = upperTriangle(X,diag=F,byrow=T)
  xld = lowerTriangle(X,diag=F,byrow=F)

  S = mean(abs(xod-xld))
  return(S=S)
}

T_comp_b3 = function(X,Y){
  xod = lowerTriangle(X,diag=F,byrow=F)
  yod = lowerTriangle(Y,diag=F,byrow=F)
  S = mean(abs(xod-yod))
  return(S=S)
}

mant_comp = function(X,Y){
  xod = upperTriangle(X,diag=F,byrow=T)
  yod = upperTriangle(Y,diag=F,byrow=T)
  C = cor(xod,yod)
  return(C=C)
}

mant_comp_b_u = function(X,Y){
  xod = upperTriangle(X,diag=F,byrow=T)
  yod = upperTriangle(Y,diag=F,byrow=T)
  C = cor(xod,yod)
  return(C=C)
}

mant_comp_b_l = function(X,Y){
  xod = lowerTriangle(X,diag=F,byrow=F)
  yod = lowerTriangle(Y,diag=F,byrow=F)
  C = cor(xod,yod)
  return(C=C)
}

mant_comp_b = function(X,Y){
  xod = upperTriangle(X,diag=T,byrow=T)
  xld = lowerTriangle(X,diag=F,byrow=T)
  x = c(xod,xld)
  yod = upperTriangle(Y,diag=T,byrow=T)
  yld = lowerTriangle(Y,diag=F,byrow=T)
  y = c(yod,yld)
  S = cor(x,y)
  return(S=S)
}


mant_comp_b2 = function(X){
  xod = upperTriangle(X,diag=F,byrow=T)
  xld = lowerTriangle(X,diag=F,byrow=F)
  S = cor(xod,xld)
  return(S=S)
}

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



RS_unif_quick = function(X,Y,nsim,sim){
  m = nrow(X)
  set.seed(4324352)
  if(sim == "norm"){
  S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  } else {
    S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  }
  R1 <- matrix(rep(1,nsim*m), nsim, m)
  R2 <- matrix(rep(1,nsim*m), nsim, m)
  for (i in 1:nrow(R1)){
    R1[i,] = (S[i,])%*%X
  }
  for (i in 1:nrow(R2)){
    R2[i,] = (S[i,])%*%Y
  }
  mvc= mean(rowSums(R1*R2)/(sqrt(rowSums(R1^2))*sqrt(rowSums(R2^2))))

  return(mvc=mvc)
}


RS_unif_magnitude_ratio = function(X,Y,nsim,sim){
  m = nrow(X)
  set.seed(4324352)
  if(sim == "norm"){
    S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  } else {
    S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  }
  
  R1 = matrix(t(apply(S,1,"%*%",X)),ncol=m,nrow=nsim)
  R2 = matrix(t(apply(S,1,"%*%",Y)),ncol=m,nrow=nsim)
  
  r1n = rowSums(R1^2)
  r2n = rowSums(R2^2)
 meanrat = sqrt(mean(r2n))/sqrt(mean(r1n))
  return(meanrat)
}


RS_unif_magnitude_sqrt = function(X,nsim,sim){
  m = nrow(X)
  set.seed(4324352)
  if(sim == "norm"){
    S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  } else {
    S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  }
  
  R1 = matrix(t(apply(S,1,"%*%",X)),ncol=m,nrow=nsim)

  
 
  r1n = rowSums(R1^2)
  meanr = mean(r1n)

  return(sqrt(meanr))
}


RS_unif_magnitude = function(X,nsim,sim){
  m = nrow(X)
  set.seed(4324352)
  if(sim == "norm"){
    S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  } else {
    S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  }
  
  R1 = matrix(t(apply(S,1,"%*%",X)),ncol=m,nrow=nsim)
  
  
  
  r1n = rowSums(R1^2)
  meanr = mean(r1n)
  return(meanr)
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

RS_m_std_beta = function(X,Y,mx,my,nsim){
  m = nrow(X)
  set.seed(493242)
  S <- matrix(runif(nsim * m, min = -2, max = 2), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  R1 <- matrix(rep(1,nsim*m), nsim, m)
  R2 <- matrix(rep(1,nsim*m), nsim, m)
  for (i in 1:nrow(R1)){
    R1[i,] = (S[i,]/mx)%*%X
  }
  for (i in 1:nrow(R2)){
    R2[i,] = (S[i,]/my)%*%Y
  }
  mvc = mean(rowMeans(R1*R2)/sqrt(rowMeans(R1^2)*rowMeans(R2^2)))
  return(mvc = mvc)
  
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
  set.seed(4604565)
  S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
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
  vcdist = rowSums(Rtop*Rbot)/(sqrt(rowSums(Rtop^2))*sqrt(rowSums(Rbot^2)))
  vectorcordist = rowMeans(Rtop*Rbot)/sqrt(rowMeans(Rtop^2)*rowMeans(Rbot^2))
  meanvectorcor = mean(rowMeans(Rtop*Rbot)/(sqrt(rowMeans(Rtop^2))*sqrt(rowMeans(Rbot^2))))
  meanevol = mean(rowMeans(R*S)/(sqrt(rowMeans(R^2))*sqrt(rowMeans(S^2))))
  evodist = rowMeans(R*S)/(sqrt(rowMeans(R^2))*sqrt(rowMeans(S^2)))
  return(list(S = S,vcd = vectorcordist,vcd2 =vcdist,mvc = meanvectorcor,r=r,corbytrait = corbytrait,Rtop = Rtop,Rbot = Rbot,evodist=evodist,meanevol=meanevol))
}










SAskewer_quick = function(X,nsim,type,sim){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  set.seed(641980)
  if(sim=="norm"){
  S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  } else {
    S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  }
  
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
  vcdist = rowSums(Rtop*Rbot)/(sqrt(rowSums(Rtop^2))*sqrt(rowSums(Rbot^2)))
  vectorcordist = rowMeans(Rtop*Rbot)/sqrt(rowMeans(Rtop^2)*rowMeans(Rbot^2))
  meanvectorcor = mean(rowSums(Rtop*Rbot)/(sqrt(rowSums(Rtop^2))*sqrt(rowSums(Rbot^2))))
  
  return(meanvectorcor)
  
}



#### sa quick sd betas
SAskewer_quick_sd = function(X,sd,nsim,type,sim){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  set.seed(641980)
  if(sim=="norm"){
    S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  } else {
    S <- matrix(runif(nsim * m, min = -1, max = 1), nsim, m)
    S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  }
  
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
    R[i,] = (0.5*(S[i,]*sd)%*%X)*sd
  }
  
  Rtop = R[,1:(nrow(X)/2)]
  Rbot = R[,((nrow(X)/2)+1):nrow(X)]
  vcdist = rowSums(Rtop*Rbot)/(sqrt(rowSums(Rtop^2))*sqrt(rowSums(Rbot^2)))
  vectorcordist = rowMeans(Rtop*Rbot)/sqrt(rowMeans(Rtop^2)*rowMeans(Rbot^2))
  meanvectorcor = mean(rowSums(Rtop*Rbot)/(sqrt(rowSums(Rtop^2))*sqrt(rowSums(Rbot^2))))
  
  return(meanvectorcor)
  
}



#### sa skewer quick compare 2 mats


SAskewer_quick_compare = function(X,Y,sdx,sdy,nsim,type){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  set.seed(641980)
  S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)

  Rx <- matrix(rep(1,nsim*m), nsim, m)
  Ry  <- matrix(rep(1,nsim*m), nsim, m)

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
  
  for (i in 1:nrow(Rx)){
    Rx[i,] = (0.5*(S[i,])%*%X)/sdx
  }
  for (i in 1:nrow(Ry)){
    Ry[i,] = (0.5*(S[i,])%*%Y)/sdy
  }
  Rtopx = Rx[,1:(nrow(X)/2)]
  Rbotx = Rx[,((nrow(X)/2)+1):nrow(X)]
  Rtopy = Ry[,1:(nrow(X)/2)]
  Rboty = Ry[,((nrow(X)/2)+1):nrow(X)]

  vcdistx = rowSums(Rtopx*Rbotx)/(sqrt(rowSums(Rtopx^2))*sqrt(rowSums(Rbotx^2)))
  vectorcordistx = rowMeans(Rtopx*Rbotx)/sqrt(rowMeans(Rtopx^2)*rowMeans(Rbotx^2))
  meanvectorcorx = mean(rowMeans(Rtopx*Rbotx)/(sqrt(rowMeans(Rtopx^2))*sqrt(rowMeans(Rbotx^2))))
  
  vcdisty = rowSums(Rtopy*Rboty)/(sqrt(rowSums(Rtopy^2))*sqrt(rowSums(Rboty^2)))
  vectorcordisty = rowMeans(Rtopy*Rboty)/sqrt(rowMeans(Rtopy^2)*rowMeans(Rboty^2))
  meanvectorcory = mean(rowMeans(Rtopy*Rboty)/(sqrt(rowMeans(Rtopy^2))*sqrt(rowMeans(Rboty^2))))

  md = meanvectorcorx-meanvectorcory
  return(c(meanvectorcorx,meanvectorcory,md))
  
}

#### chen and Houle selection response in antagonistic and concordant space
AC_skewers_full = function(X,nsim){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  
  I = diag((m/2))
  Q = 0.5*(rbind(cbind(I,I),cbind(I,-I)))
  Gca = Q%*%X%*%Q
  
  
  S0 = matrix(rep(0,nsim*(m/2)),nsim,(m/2))
  S <- matrix(runif(nsim * (m/2), min = -2, max = 2), nsim, (m/2))
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, (m/2))
  Sc = cbind(S,S0)
  Sa = cbind(S0,S)
  
  R <- matrix(rep(1,nsim*m), nsim, m)
  Rc = R
  Ra = R
  for (i in 1:nrow(Rc)){
    Rc[i,] = (Sc[i,])%*%Gca
  }
  for (i in 1:nrow(Ra)){
    Ra[i,] = (Sa[i,])%*%Gca
  }
  Rcc = Rc[,1:(nrow(X)/2)]
  Rca = Rc[,((nrow(X)/2)+1):nrow(X)]
  Rac = Ra[,1:(nrow(X)/2)]
  Raa = Ra[,((nrow(X)/2)+1):nrow(X)]
  
  return(list(Rcc = Rcc,Rca=Rca,Rac=Rac,Raa=Raa,Sa=Sa,Sc=Sc,Ra=Ra,Rc=Rc))
}


AC_skewers_quick = function(X,nsim){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  
  I = diag((m/2))
  Q = 0.5*(rbind(cbind(I,I),cbind(I,-I)))
  Gca = Q%*%X%*%Q
  
  
  S0 = matrix(rep(0,nsim*(m/2)),nsim,(m/2))
  S <- matrix(runif(nsim * (m/2), min = -2, max = 2), nsim, (m/2))
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, (m/2))
  Sc = cbind(S,S0)
  Sa = cbind(S0,S)
  
  R <- matrix(rep(1,nsim*m), nsim, m)
  Rc = R
  Ra = R
  for (i in 1:nrow(Rc)){
    Rc[i,] = (Sc[i,])%*%Gca
  }
  for (i in 1:nrow(Ra)){
    Ra[i,] = (Sa[i,])%*%Gca
  }
  Rcc = Rc[,1:(nrow(X)/2)]
  Rca = Rc[,((nrow(X)/2)+1):nrow(X)]
  Rac = Ra[,1:(nrow(X)/2)]
  Raa = Ra[,((nrow(X)/2)+1):nrow(X)]
  
  m_mag_aa = mean(sqrt(rowSums(Raa^2)))
  m_mag_ac = mean(sqrt(rowSums(Rac^2)))
  m_mag_ca = mean(sqrt(rowSums(Rca^2)))
  m_mag_cc = mean(sqrt(rowSums(Rcc^2)))
  return(c(m_mag_aa=m_mag_aa,m_mag_ac=m_mag_ac,m_mag_cc=m_mag_cc,m_mag_ca=m_mag_ca))
}


#### measure B matrix symmetry
B_sym = function(X){
  m = nrow(X)
  B = X[1:(0.5*nrow(X)),((0.5*nrow(X))+1):nrow(X)]
  r = (sum(diag(t(t(B))%*%B)))/(norm(B,type="F")*norm(t(B),type="F"))
  return(r=r)
}
Sym = function(X,Y){
  m = nrow(X)
  r = (sum(diag(t(X)%*%Y)))/(norm(X,type="F")*norm(Y,type="F"))
  return(r=r)
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


