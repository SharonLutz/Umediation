generate_data_matrix <- function(
  n=100,interact=FALSE,
  Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",
  muC=0,muU=0,varA=1,varM=1,varY=1,varC=1,varU=1,
  gamma0=0,gammaC=0,gammaU=0,
  alpha0=0,alphaA=0,alphaC=0,alphaU=0,
  beta0=0,betaA=0,betaM=0,betaI=0,betaC=0,betaU=0,
  nSim=250,seed=1,nBoot=400){
  
  data_matrix = list()
  matC = vector(mode="numeric")
  for(si in 1:nSim){ #loop through all simulations
    set.seed(seed+(si-1)) #set the seed to ensure reproducibility
    
    # simStep<-1
    # if(floor(si/simStep)==ceiling(si/simStep)){print(paste("simulation",si, "of",nSim))}
    
    #######################################
    # Generate the unmeasured confounder,
    # measured confounder C, exposure A,
    # mediator M, and outcome Y
    #######################################
    
    # generate unmeasured confounder U
    U<-matrix(0,nrow=n,ncol=length(Utype))
    for(iu in 1:length(Utype)){
      if(Utype[iu]=="C"){
        if(varU[iu]==0|varU[iu]<0){stop(paste("Error: Variance of U",iu," (varU) is 0 or negative.",sep=""))}
        U[,iu]<-rnorm(n,muU[iu],sqrt(varU[iu]))
      }
      if(Utype[iu]=="D"){
        if(muU[iu]<0|muU[iu]==0|muU[iu]==1|muU[iu]>1){stop(paste("Error: Prob(U",iu,"=1) (i.e. muU",iu,") is less than or equal to zero or greater than or equal to 1.",sep=""))}
        U[,iu]<-rbinom(n,1,muU[iu])
      }
      colnames(U)<-paste("U",length(Utype),Utype,sep="")
    }
    
    # generate measured confounder C
    CC<-matrix(0,nrow=n,ncol=length(Ctype))
    for(ic in 1:length(Ctype)){
      if(Ctype[ic]=="C"){
        if(varC[ic]==0|varC[ic]<0){stop(paste("Error: Variance of C",ic," (varC) is 0 or negative.",sep=""))}
        CC[,ic]<-rnorm(n,muC[ic],sqrt(varC[ic]))
      }
      if(Ctype[ic]=="D"){
        if(muC[ic]<0|muC[ic]==0|muC[ic]==1|muC[ic]>1){stop(paste("Error: Prob(C",ic,"=1) is less than or equal to zero or greater than or equal to 1.",sep=""))}
        CC[,ic]<-rbinom(n,1,muC[ic])
      }
      colnames(CC)<-paste("C",length(Ctype),Ctype,sep="")
    }
    
    # generate exposure A as a function of C and U
    # logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C
    muA<-gamma0+CC%*%gammaC+U%*%gammaU
    if(Atype=="C"){
      A<-rnorm(n,muA,sqrt(varA)) # normally distributed mediator
      if(varA==0|varA<0){stop("Error: Variance of M (varM) is 0 or negative.")}
    }
    if(Atype=="D"){
      Ap<-exp(muA)/(1+exp(muA))
      if(max(Ap)=="NaN"|is.na(max(Ap))){stop("Error: Prob(A=1) is too close to 0 or 1. This can occur if the absolute value of gamma0 is too large.")}
      A<-rbinom(n,1,Ap)# binary mediator
      if(length(A[A==1])/length(A)<(2/n)|length(A[A==1])/length(A)>(1-2/n)){
        stop("Error: there was no enough varaibality in A (i.e. Prob(A=1)=1 or Prob(A=1)=0). Considering increasing the sample size n or using mean centered confounders (i.e. muC=0 and muU=0).")
      }
    }
    
    # generate mediator M
    # Logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U)
    muM<-alpha0+alphaA*A+CC%*%alphaC+U%*%alphaU
    if(Mtype=="C"){
      M<-rnorm(n,muM,sqrt(varM)) # normally distributed mediator
      if(varM==0|varM<0){stop("Error: Variance of M (varM) is 0 or negative.")}
    }
    if(Mtype=="D"){
      Mp<-exp(muM)/(1+exp(muM))
      if(max(Mp)=="NaN"|is.na(max(Mp))){stop("Error: Prob(M=1) is too close to 0 or 1. This can occur if the absolute value of alpha0 is too large.")}
      M<-rbinom(n,1,Mp)# binary mediator
      if(length(M[M==1])/length(M)<(5/n)|length(M[M==1])/length(M)>(1-5/n)){
        stop("Error: there was no enough varaibality in M (i.e. Prob(M=1)=1 or Prob(M=1)=0). Considering increasing the sample size n or using mean centered confounders (i.e. muC=0 and muU=0).")
      }
    }
    
    
    # generate outcome Y
    # logit(P(Y=1)) or E[Y]=beta0+betaA*A+betaM*M+betaC*C+betaU*U
    if(interact==TRUE){ muY<-beta0+betaA*A+betaM*M+betaI*A*M+CC%*%betaC+U%*%betaU}
    if(interact==FALSE){ muY<-beta0+betaA*A+betaM*M+CC%*%betaC+U%*%betaU}
    if(Ytype=="C"){
      if(varY==0|varY<0){stop("Error: Variance of Y (varY) is 0 or negative.")}
      Y<-rnorm(n,muY,sqrt(varY)) # normally distributed outcome
    }
    if(Ytype=="D"){
      Yp<-exp(muY)/(1+exp(muY))
      if(max(Yp)=="NaN"|is.na(max(Yp))){stop("Error: Prob(Y=1) is too close to 0 or 1. This can occur if the absolute value of beta0 is too large.")}
      Y<-rbinom(n,1,Yp) # binary outcome
      if(length(Y[Y==1])/length(Y)<(5/n)|length(Y[Y==1])/length(Y)>(1-5/n)){
        stop("Error: there was no enough varaibality in Y (i.e. Prob(Y=1)=1 or Prob(Y=1)=0). Considering increasing the sample size n or using mean centered confounders (i.e. muC=0 and muU=0).")
      }
    }
    
    #######################################
    # Check for collinearity
    #######################################
    
    if(length(Utype)>1|length(Ctype)>1){
      if(max(car::vif(lm(Y~A+M+CC+U))[,3])>10 & interact==FALSE ){stop("Error: There is strong evidence of collinearity (i.e. car::vif>10) when A,M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(car::vif(lm(Y~A+M+A*M+CC+U))[,3])>10 & interact==TRUE ){stop("Error: There is strong evidence of collinearity (i.e. car::vif>10) when A,M,A*M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(car::vif(lm(M~A+CC+U))[,3])>10){stop("Error: There is strong evidence of collinearity (i.e. car::vif>10) when A,C,U are regressed on M. Consider reducing the absolute value of large gammas or alphas.")}
    }
    
    if(length(Utype)==1&length(Ctype)==1){
      if(max(car::vif(lm(Y~A+M+CC+U)))>10 & interact==FALSE ){stop("Error: There is strong evidence of collinearity (i.e. car::vif>10) when A,M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(car::vif(lm(Y~A+M+A*M+CC+U)))>10 & interact==TRUE ){stop("Error: There is strong evidence of collinearity (i.e. car::vif>10) when A,M,A*M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
      
      if(max(car::vif(lm(M~A+CC+U)))>10){stop("Error: There is strong evidence of collinearity (i.e. car::vif>10) when A,C,U are regressed on M. Consider reducing the absolute value of large gammas or alphas.")}
    }
    
    #######################################
    # Correlation for U,C,A,M,Y
    #######################################
    
    matA<-cbind(A,M,Y,CC,U)
    if(length(matC) > 0){
      matC <- matC + cor(matA)
    } else {
      matC <- cor(matA)
    }
    colnames(matC)<-c("A","M","Y",paste("C",c(1:length(Ctype)),sep=""),paste("U",c(1:length(Utype)),sep=""))
    rownames(matC)<-c("A","M","Y",paste("C",c(1:length(Ctype)),sep=""),paste("U",c(1:length(Utype)),sep=""))
    
    # TODO: store the input data into output matrix
    
    data_matrix[[si]] = list(
      interact=interact,
      Mtype=Mtype,
      Ytype=Ytype,
      M=M,
      Y=Y,
      A=A,
      CC=CC,
      U=U,
      nBoot=nBoot,
      RandStateParam = .Random.seed
    )
    
  } #end of simulation loop
  return(list(
    matC=(round((matC / nSim), digits=2)),
    data=data_matrix
  ))
}