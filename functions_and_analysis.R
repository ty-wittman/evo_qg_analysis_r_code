
### script load in all functions#####
## modified t method comparison
T_comp = function(X,Y){
  xod = upperTriangle(X,diag=F,byrow=T)
  yod = upperTriangle(Y,diag=F,byrow=T)
  S = abs(xod-yod)
  return(S=S)
}

#### mantel correlation

mant_comp = function(X,Y){
  xod = upperTriangle(X,diag=F,byrow=T)
  yod = upperTriangle(Y,diag=F,byrow=T)
  C = cor(xod,yod)
  return(C=C)
}


### function to do random skewers analysis. X and Y are symmetric covariance matrices, nsim is the number of simulated selection skewers

RS = function(X,Y,nsim){
  m = nrow(X)
  S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
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

### function to calculate average between sex genetic correlation, can be used in conjunction with the sampling distribution of the full G matrices to calculate a sampling distribution of mean between sex genetic correlation. 
### this is not in a generalized form, X must be a 10x10 genetic correlation matrix with the upper right and lower left 25 cells containing the between sex components.
rmf_m = function(X){
  S = mean(diag(X[6:10,1:5]))
  return(S)
}

## this function calculates the average of the differences between paired between sex genetic correlations for two matricies
### can be used to generate a sampling distribution of mean differences
### this is not in a generalized form, X and Y must be 10x10 genetic correlation matrices with the upper right and lower left 25 cells containing the between sex components.

rmf_d = function(X,Y){
  S = mean(diag(X[6:10,1:5])-diag(Y[6:10,1:5]))
  return(S)
}

## This function performs sexually antagonistic skewers analysis.  X is the full G matrix which includes B, nsim is the number of random sexually antagonistic skewers. 
### not generalizable at the moment, in making the skewers sexually antagonistic the function assumes a 10x10 matrix.
SAskewers = function(X,nsim){
  m = nrow(X)
  if(!all(sapply(dim(X),"==",m)))
    stop("X should be a square matrix")
  if((nrow(X)%%2)!=0)
    stop("X should be divisable by two")
  S <- matrix(rnorm(nsim * m, mean = 0, sd = 1), nsim, m)
  S <- S/matrix(sqrt(rowSums(S^2)), nsim, m)
  for (i in seq(1,5)){
        S[,i+5]= ifelse(S[,i] >0 & S[,i+5] > 0 , -S[,i+5],S[,i+5] )
        S[,i+5] = ifelse(S[,i] <0 & S[,i+5] <0 , -S[,i+5],S[,i+5] )
     }
  R <- matrix(rep(1,nsim*m), nsim, m)
  for (i in 1:nrow(R)){
    R[i,] = 0.5*(S[i,]%*%X)
  }
  
  Rtop = R[,1:(ncol(R)/2)]
  Rbot = R[,((ncol(R)/2)+1):ncol(R)]
  meanvectorcor = mean(rowSums(Rtop*Rbot)/(sqrt(rowSums(Rtop^2))*sqrt(rowSums(Rbot^2))))
  
  return(meanvectorcor)
  
}













##### analysis code ######

### packages used
library(DescTools)
library(cowplot)
library(forcats)
library(plyr)
library(dplyr)
library(tidyr)
library(GGally)
library(gtable)
library(stringr)
library(plyr)
library(tidyr)
library(gdata)
library(snow)
library(doParallel)
library(foreach)
library(snow)
library(doParallel)
library(foreach)
library(Matrix)
library(Rfast)
library(tidyr)
library(matrixcalc)
library(LaplacesDemon)

#### data is avaliable on dryat at https://doi.org/10.5061/dryad.1rn8pk0tf
### read in phenotypic data

# control males
cmpmat = read.csv("path to file/control_male_phenotypic_matrix.csv")
cmpmat = as.matrix(cmpmat[,2:6])
cmsd = sqrt(diag(cmpmat))

#control females
cfpmat = read.csv("path to file/control_female_phenotypic_matrix.csv")
cfpmat = as.matrix(cfpmat[,2:6])
cfsd = sqrt(diag(cfpmat))

#testosterone males
tmpmat = read.csv("path to file/testosterone_male_phenotypic_matrix.csv")
tmpmat = as.matrix(tmpmat[,2:6])
tmsd = sqrt(diag(tmpmat))
#testosterone females
tmpmat = read.csv("path to file/testosterone_female_phenotypic_matrix.csv")
tmpmat = as.matrix(tmpmat[,2:6])
tmsd = sqrt(diag(tmpmat))

### these will be used for analysis and to standardize the genetic variance-covariance matrices
### both unstandardized and standardized G matrices are avaliable on dryad
### read in the genetic variance-covariance matrices

#control male
cmgmat = read.csv("path to file/Control_male_distribution_unstandardized.csv")
cmgmat = as.matrix(cmgmat[,2:6])

#control female
cfgmat = read.csv("path to file/Control_female_distribution_unstandardized.csv")
cfgmat = as.matrix(cfgmat[,2:6])

#testosterone male
tmgmat = read.csv("path to file/Testosterone_male_distribution_unstandardized.csv")
tmgmat = as.matrix(tmgmat[,2:6])

#testosterone female
tmgmat = read.csv("path to file/Testosterone_female_distribution_unstandardized.csv")
tmgmat = as.matrix(tmgmat[,2:6])



#### standardize the G matrices by the phenotypic variances

vs_cmmat = cmgmat/(sdbm%*%t(sdbm))
vs_cfmat = cfmgat/(sdbf%*%t(sdbf))
vs_tmmat = tmgmat/(sdtm%*%t(sdtm))
vs_tfmat = tfgmat/(sdtf%*%t(sdtf))



### do all pairwise comparisons between variance standardized best estimate G matrices


std_cftf = RS_quick(vs_cfmat,vs_tfmat,10000)
std_cfcm = RS_quick(vs_cfmat,vs_cmmat,10000)
std_cftm = RS_quick(vs_cfmat,vs_tmmat,10000)
std_tfcm = RS_quick(vs_tfmat,vs_cmmat,10000)
std_tftm = RS_quick(vs_tfmat,vs_tmmat,10000)
std_cmtm = RS_quick(vs_cmmat,vs_tmmat,10000)


v_std_g_comp = data.frame(std_cftf,std_cfcm,std_cftm,std_tfcm,std_tftm,std_cmtm)
v_std_g_comp




### read in sampling distribution for each sex and treatment group. 
### these are in CSV files that contain a total of 10,000 matrices, the matrices are fully filled in, so the file needs to be sliced every 10 rows. 
### lot of ways you can do this. I like to use arrays 
### in the data file there is a column giving the number of each matrix, 1 through 10000, this may be useful for splitting the file if you prefor to use a different format
samp_dist_cm = read.csv("path to file/Control_male_distribution_unstandardized.csv")
samp_dist_cm = as.matrix(samp_dist_cm[,2:6])
samp_dist_cm = array(aperm(samp_dist_cm),dim=c(5,5,10000))

samp_dist_cf = read.csv("path to file/Control_female_distribution_unstandardized.csv")
samp_dist_cf = as.matrix(samp_dist_cf[,2:6])
samp_dist_cf = array(aperm(samp_dist_cf),dim=c(5,5,10000))

samp_dist_tm = read.csv("path to file/Testosterone_male_distribution_unstandardized.csv")
samp_dist_tm = as.matrix(samp_dist_tm[,2:6])
samp_dist_tm = array(aperm(samp_dist_tm),dim=c(5,5,10000))

samp_dist_tf = read.csv("path to file/Testosterone_female_distribution_unstandardized.csv")
samp_dist_tf = as.matrix(samp_dist_tm[,2:6])
samp_dist_tf = array(aperm(samp_dist_tm),dim=c(5,5,10000))




### now we need to standardize each matrix in the sampling distribution. 


samp_dist_v_std_cm = array(apply(samp_dist_cm,3,function(x) round((x/(sdcm%*%t(sdcm))),10)),dim=c(5,5,10000))

samp_dist_v_std_cf = array(apply(samp_dist_cf,3,function(x) round((x/(sdcf%*%t(sdcf))),10)),dim=c(5,5,10000))
                                  
samp_dist_v_std_tm = array(apply(samp_dist_tm,3,function(x) round((x/(sdtm%*%t(sdtm))),10)),dim=c(5,5,10000))
                                 
samp_dist_v_std_tf = array(apply(samp_dist_tf,3,function(x) round((x/(sdtf%*%t(sdtf))),10)),dim=c(5,5,10000))

                                 
### now to generate the null distribution for each matrix by comparing the best estimate to each matrix in its sampling distribution with random skewers, mantel comparison and modified T method
                                 

### Running this parallel saves a lot of time. 
## i use snow, doparrallel and foreach, 
cl = makeCluster(4)
registerDoParallel(cl)
### need to read your best estimate matrices and functions into the cluster
clusterExport(cl, list("vs_cmmat","vs_cfmat","vs_tfmat","vs_tmmat","RS"))
                                 
### use parralell apply function to get null sampling distribution, using 10000 random skewers
vs_std_cm=parApply(cl,v_std_cmsg,3,function(x) RS(x,vs_cmmat,nsim=10000))                        
                                 
vs_std_cf=parApply(cl,v_std_cfsg,3,function(x) RS(x,vs_cfmat,nsim=10000))                                  

vs_std_tm=parApply(cl,v_std_tmsg,3,function(x) RS(x,vs_tmmat,nsim=10000)) 

vs_std_tf=parApply(cl,v_std_tfsg,3,function(x) RS(x,vs_tfmat,nsim=10000))     
                   
                   
## We have all our samping distributions, now we can get information about the distribution, The mode, the lower 5% cutoff, and the count below our best estimate

cm_mode =  LaplacesDemon::Mode(vs_std_cm)     
low_cut_off_cm = quantile(vs_std_cm,0.05)
                   
plyr::count(vs_std_cm<v_std_g_comp["std_cfcm"])
plyr::count(vs_std_cm<v_std_g_comp["std_tfcm"])
plyr::count(vs_std_cm<v_std_g_comp["std_cmtm"])

                   
cf_mode =  LaplacesDemon::Mode(vs_std_cf)     
low_cut_off_cf = quantile(vs_std_cf,0.05)
                   
plyr::count(vs_std_cf<v_std_g_comp["std_cftf"])
plyr::count(vs_std_cf<v_std_g_comp["std_cfcm"])
plyr::count(vs_std_cf<v_std_g_comp["std_cftm"])
                   
                   
tm_mode =  LaplacesDemon::Mode(vs_std_tm)     
low_cut_off_tm = quantile(vs_std_tm,0.05)
                   
plyr::count(vs_std_tm<v_std_g_comp["std_cftm"])
plyr::count(vs_std_tm<v_std_g_comp["std_tftm"])
plyr::count(vs_std_tm<v_std_g_comp["std_cmtm"])
                        
tf_mode =  LaplacesDemon::Mode(vs_std_tf)     
low_cut_off_tf = quantile(vs_std_tf,0.05)
                   
plyr::count(vs_std_tf<v_std_g_comp["std_cftf"])
plyr::count(vs_std_tf<v_std_g_comp["std_tfcm"])
plyr::count(vs_std_tf<v_std_g_comp["std_tftm"])
    

####  To do analysis on unstandardized sampling distribution just replace the standardized one with the unstandardized one

                   
 ### for t method and mantel correlation we need to convert the genetic covariance matrices to genetic correlation matrices
 ## the t method and mantel corerlation we use just compares the off-diagonals so measures matrix shape not matrix size
 cor_cmgmat = cov2cor(cmgmat)
 cor_cfgmat = cov2cor(cfgmat)
 cor_tmgmat = cov2cor(tmgmat)                  
 cor_tfgmat = cov2cor(tfgmat)
                   
  ### convert sampling distributions from covariance to correlation
samp_dist_cor_cm = array(apply(samp_dist_cm,3,function(x) cov2cor(x),dim=c(5,5,10000))

samp_dist_cor_cf = array(apply(samp_dist_cf,3,function(x) cov2cor(x),dim=c(5,5,10000))
                                  
samp_dist_cor_tm = array(apply(samp_dist_tm,3,function(x) cov2cor(x),dim=c(5,5,10000))
                                 
samp_dist_cor_tf = array(apply(samp_dist_tf,3,function(x) cov2cor(x),dim=c(5,5,10000))
                   
                               
## now to do the t method and mantel cor comparison just take the code from the RS section and replace RS with T_comp or mant_comp
                               
                   
#### now for sexually antagonistic skewers.                    
### read in matrices used to standardize the full G with B matrix, this standardization matrix is the cross product of the phenotypic standard deviations       
                   
cmcfsd = read.csv("path to file/control_male_control_female_cross_product_standard_deviation_for_variance_standardization.csv")
cmcfsd = as.matrix(cmcfsd)


cmtfsd= read.csv("path to file/control_male_testosterone_female_cross_product_standard_deviation_for_variance_standardization.csv")
cmtfsd = as.matrix(cmtfsd)
  
tmtfsd= read.csv("path to file/testosterone_male_testosterone_female_cross_product_standard_deviation_for_variance_standardization.csv")
tmtfsd = as.matrix(cmtfsd)
                   

### read in genetic variance-covariance matrices with B estimated for two groups
### standardized as well as unstandardized matrices are given on dryad.
                   
cmcfgmat = read.csv("path to file/control_male_control_female_full_G_unstandardized.csv")
cmcfgmat = as.matrix(cmcfgmat)   
                   
cmtfgmat = read.csv("path to file/control_male_testosterone_female_full_G_unstandardized.csv")
cmtfgmat = as.matrix(cmtfgmat)               
                   
 tmtfgmat = read.csv("path to file/testosterone_male_testosterone_female_full_G_unstandardized.csv")
tmtfgmat = as.matrix(tmtfgmat)   
                   
                   
### standardize G with B using phenotypic variances
                   
                   
vs_cmcfmat = cmcfgmat/cmcfsd                  
 
vs_cmtfmat = cmtfgmat/cmtfsd  
      
vs_tmtfmat = tmtfgmat/tmtfsd                                 
                   
#### measure SA response for each G with B matrix
                   
                   
SA_cmcf = SAskewers(vs_cmcfmat,10000)
                   
SA_cmtf = SAskewers(vs_cmtfmat,10000)        
                        
SA_tmtf = SAskewers(vs_tmtfmat,10000)                                
                   
### read in sampling distribution of each G with B
                   
samp_dist_cmcf = read.csv("path to file/control_male_control_female_distribution_unstandardized.csv")
samp_dist_cmcf = as.matrix(samp_dist_cmcf[,2:11])
samp_dist_cmcf = array(aperm(samp_dist_cmcf),dim=c(10,10,10000))
                   
samp_dist_cmtf = read.csv("path to file/control_male_testosterone_female_distribution_unstandardized.csv")
samp_dist_cmtf = as.matrix(samp_dist_cmtf[,2:11])
samp_dist_cmtf = array(aperm(samp_dist_cmcf),dim=c(10,10,10000))    
                   
                   
samp_dist_tmtf = read.csv("path to file/testosterone_male_testosterone_female_distribution_unstandardized.csv")
samp_dist_tmtf = as.matrix(samp_dist_tmtf[,2:11])
samp_dist_tmtf = array(aperm(samp_dist_tmcf),dim=c(10,10,10000))                     
                   
                   
                   
                   
### standardize matrices in sampling distribution
                   
samp_dist_v_std_cmcf = array(apply(samp_dist_cmcf,3,function(x) round(x/cmcfsd,10),dim=c(10,10,10000))

samp_dist_v_std_cmtf = array(apply(samp_dist_cmtf,3,function(x) round(x/cmtfsd,10),dim=c(10,10,10000))        

samp_dist_v_std_cmtf = array(apply(samp_dist_cmtf,3,function(x) round(x/cmtfsd,10),dim=c(10,10,10000))        

                                   
## set up cluster                        
cl = makeCluster(#number of clusters you want to use)
registerDoParallel(cl)
### need to read your best estimate matrices and functions into the cluster
clusterExport(cl, "vs_cmcfmat","vs_cmtfmat","vs_tmtfmat","SAskewers")                    
                                   
                                   
##### get SAskewers estimate for each matrix in sampling distribution, create null distribution                                   
                                   
vs_std_cmcf=parApply(cl,samp_dist_v_std_cmcf,3,function(x) SAskewers(x,nsim=10000))                        
                                 
vs_std_cmtf=parApply(cl,samp_dist_v_std_cmtf,3,function(x) SAskewers(x,nsim=10000))                                    
 
vs_std_tmtf=parApply(cl,samp_dist_v_std_tmtf,3,function(x) SAskewers(x,nsim=10000))                     
                                   
cmcf_mode =  LaplacesDemon::Mode(vs_std_cmcf)     
low_cut_off_cmcf = quantile(vs_std_cmcf,0.95)
### is correlation between control male and testosterone female responses to sexually antagonistic selection significantly greater than that of control males and control females
plyr::count(vs_std_cmcf>SA_cmtf)
### is correlation between testosterone male and testosterone female responses to sexually antagonistic selection significantly greater than that of control males and control females
plyr::count(vs_std_cmcf>SA_tmtf)

                                   
cmtf_mode =  LaplacesDemon::Mode(vs_std_cmtf)     
low_cut_off_cmtf = quantile(vs_std_cmtf,0.05)
### is the correlation between control male and control female responses to sexually antagonistic selection significantly lower than that of control males and testosterone females
                     
plyr::count(vs_std_cmtf<SA_cmcf)           

  
                     
cmtf_mode =  LaplacesDemon::Mode(vs_std_tmtf)     
low_cut_off_tmtf = quantile(vs_std_tmtf,0.05)
### is the correlation between control male and control female responses to sexually antagonistic selection significantly lower than that of testosterone males and testosterone females
                     
plyr::count(vs_std_tmtf<SA_cmcf)   

                     
### test if the average strength of between sex genetic correlations differ between control males and control females, control males and testosterone females,  and testosterone males and testosterone females. 
                     
                     
                     
                     
