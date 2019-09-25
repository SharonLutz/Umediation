
Umediation<-function(n=100,Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",
                     interact=FALSE,muC=0,varC=1,muU=0,varU=1,
                     gamma0=0,gammaC=0,gammaU=0,varA=1,
                     alpha0=0,alphaA=0,alphaC=0,alphaU=0,varM=1,
                     beta0=0,betaA=0,betaM=0,betaI=0,betaC=0,betaU=0,varY=1,
                     alpha=0.05,nSim=250,nBoot=400,seed=1,atreat=1,acontrol=0){
    
    library(mediation) #load the mediation library
    library(car) #load the car library
    
    #######################################
    # Check Input for Errors
    #######################################
    ErrorCheck(n,Atype,Mtype,Ytype,Ctype,Utype,interact,muC,varC,muU,varU,gamma0,gammaC,gammaU,varA,alpha0,alphaA,alphaC,alphaU,varM,beta0,betaA,betaM,betaI,betaC,betaU,varY,alpha,nSim,nBoot,seed,atreat,acontrol)
   
    #######################################
    # Create Results Matrix
    #######################################
   
        Results<-matrix(0,nrow=12,ncol=1)
        rownames(Results)<-c("Prop. of simulations w/ significant ACME excluding U","Prop. of simulations w/ significant ACME including U","Prop. of simulations where conclusions based on ACME match","Average ACME excluding U","Average ACME including U","Average absolute difference of ACME including U minus ACME excluding U","Prop. of simulations w/ significant ADE excluding U","Prop. of simulations w/ significant ADE including U","Prop. of simulations where conclusions based on ADE match","Average ADE excluding U","Average ADE including U","Average absolute difference of ADE including U minus ADE excluding U")
        
    #correlation 
    corA<-matrix(0,nrow=(3+ncol(U)+ncol(CC)),ncol=(3+ncol(U)+ncol(CC)))
    
    set.seed(seed)
    #######################################
    # Loop through all simulations
    #######################################
    
    for(si in 1:nSim){ #loop through all simulations
        #set.seed(seed+(si-1)) #set the seed to ensure reproducibility
        
        simStep<-5
        if(floor(si/simStep)==ceiling(si/simStep)){print(paste("simulation",si, "of",nSim))}
        
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
        if(max(vif(lm(Y~A+M+CC+U))[,3])>10 & interact==FALSE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
        
        if(max(vif(lm(Y~A+M+A*M+CC+U))[,3])>10 & interact==TRUE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,A*M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
        
        if(max(vif(lm(M~A+CC+U))[,3])>10){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,C,U are regressed on M. Consider reducing the absolute value of large gammas or alphas.")}
         }
       
       if(length(Utype)==1&length(Ctype)==1){
        if(max(vif(lm(Y~A+M+CC+U)))>10 & interact==FALSE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
        
        if(max(vif(lm(Y~A+M+A*M+CC+U)))>10 & interact==TRUE ){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,M,A*M,C,U are regressed on Y. Consider reducing the absolute value of large gammas, alphas, or betas.")}
        
        if(max(vif(lm(M~A+CC+U)))>10){stop("Error: There is strong evidence of collinearity (i.e. VIF>10) when A,C,U are regressed on M. Consider reducing the absolute value of large gammas or alphas.")}
        }
        
        #######################################
        # Correlation for U,C,A,M,Y
        #######################################
        
        matA<-cbind(A,M,Y,CC,U)
        corA<-corA+cor(matA)
      print(cor(matA))
               
        #######################################
        # Mediation analysis w/ and w/out U
        #######################################
        
        # models without U
        if(Mtype=="C"){med.fit <-lm(M~A+CC)} #fit the model for the mediator
        if(Mtype=="D"){med.fit <-glm(M~A+CC,family=binomial(link = "logit"))}
        if(interact==FALSE){
            if(Ytype=="C"){out.fit<-lm(Y~M+A+CC)} #fit the model for the outcome
            if(Ytype=="D"){out.fit<-glm(Y~M+A+CC,family=binomial(link = "logit"))}
        }
        if(interact==TRUE){
            if(Ytype=="C"){out.fit<-lm(Y~M+A+CC+A*M)} #fit the model for the outcome
            if(Ytype=="D"){out.fit<-glm(Y~M+A+CC+A*M,family=binomial(link = "logit"))}
        }
        med.out <- mediate(med.fit, out.fit, treat = "A", mediator = "M",sim=nBoot)
        
        # models with U
        if(Mtype=="C"){med.fitU <-lm(M~A+CC+U)} #fit the model for the mediator
        if(Mtype=="D"){med.fitU <-glm(M~A+CC+U,family=binomial(link = "logit"))}
        if(interact==FALSE){
            if(Ytype=="C"){out.fitU<-lm(Y~M+A+CC+U)} #fit the model for the outcome
            if(Ytype=="D"){out.fitU<-glm(Y~M+A+CC+U,family=binomial(link = "logit"))}
        }
        if(interact==TRUE){
            if(Ytype=="C"){out.fitU<-lm(Y~M+A+CC+U+A*M)} #fit the model for the outcome
            if(Ytype=="D"){out.fitU<-glm(Y~M+A+CC+U+A*M,family=binomial(link = "logit"))}
        }
        med.outU <- mediate(med.fitU, out.fitU, treat = "A", mediator = "M",sim=nBoot)
        
        #Results
            Results["Average ACME excluding U",1]<-Results["Average ACME excluding U",1]+summary(med.out)$d.avg
            Results["Average ACME including U",1]<-Results["Average ACME including U",1]+summary(med.outU)$d.avg
            Results["Average absolute difference of ACME including U minus ACME excluding U",1]<-Results["Average absolute difference of ACME including U minus ACME excluding U",1]+abs(summary(med.out)$d.avg-summary(med.outU)$d.avg)
            
            Results["Average ADE excluding U",1]<-Results["Average ADE excluding U",1]+summary(med.out)$z.avg
            Results["Average ADE including U",1]<-Results["Average ADE including U",1]+summary(med.outU)$z.avg
            Results["Average absolute difference of ADE including U minus ADE excluding U",1]<-Results["Average absolute difference of ADE including U minus ADE excluding U",1]+abs(summary(med.out)$z.avg-summary(med.outU)$z.avg)
            
            if(summary(med.out)$d.avg.p<alpha){Results["Prop. of simulations w/ significant ACME excluding U",1]<- Results["Prop. of simulations w/ significant ACME excluding U",1]+1}
            if(summary(med.outU)$d.avg.p<alpha){Results["Prop. of simulations w/ significant ACME including U",1]<- Results["Prop. of simulations w/ significant ACME including U",1]+1}
            if(((summary(med.out)$d.avg.p<alpha)&(summary(med.outU)$d.avg.p<alpha))|((summary(med.out)$d.avg.p>alpha)&(summary(med.outU)$d.avg.p>alpha))){Results["Prop. of simulations where conclusions based on ACME match",1]<- Results["Prop. of simulations where conclusions based on ACME match",1]+1}
            
            if(summary(med.out)$z.avg.p<alpha){Results["Prop. of simulations w/ significant ADE excluding U",1]<- Results["Prop. of simulations w/ significant ADE excluding U",1]+1}
            if(summary(med.outU)$z.avg.p<alpha){Results["Prop. of simulations w/ significant ADE including U",1]<- Results["Prop. of simulations w/ significant ADE including U",1]+1}
            if(((summary(med.out)$z.avg.p<alpha)&(summary(med.outU)$z.avg.p<alpha))|((summary(med.out)$z.avg.p>alpha)&(summary(med.outU)$z.avg.p>alpha))){Results["Prop. of simulations where conclusions based on ADE match",1]<- Results["Prop. of simulations where conclusions based on ADE match",1]+1}
            
        
    } #end of simulation loop
    
    #######################################
    # Results
    #######################################
    matC<-round(corA/nSim,digits=2)
    colnames(matC)<-c("A","M","Y",paste("C",c(1:length(Ctype)),sep=""),paste("U",c(1:length(Utype)),sep=""))
    rownames(matC)<-c("A","M","Y",paste("C",c(1:length(Ctype)),sep=""),paste("U",c(1:length(Utype)),sep=""))
    
    matR<-Results/nSim
    listA<-list(matR,matC,paste("Warning: correlations are only valid if at least one of the variables is normally distributed."))
    names(listA)<-c("Results","Correlations_Between_Variables","Warning")
    listA
    
} #end of function
