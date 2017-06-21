
ErrorCheck<-function(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",interact=FALSE,muC=0,varC=1,muU=0,varU=1,gamma0=0,gammaC=0,gammaU=0,varA=1,alpha0=0,alphaA=0,alphaC=0,alphaU=0,varM=1,beta0=0,betaA=0,betaM=0,betaI=0,betaC=0,betaU=0,varY=1,alpha=0.05,nSim=300,nBoot=500,seed=1,atreat=1,acontrol=0){
    
    ##################
    # Models of exposure A,
    # mediator M, outcome Y
    ##################
    # logit(P(A=1)) or E[A]=gamma0+gammaC*C+gammaU*U
    # Logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U)
    # logit(P(Y=1)) or E[Y]=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U

    ##################
    # check A,M,Y,C,U type vectors are only C or D
    ##################
    if(Atype!="C"& Atype!="D"){stop("Error: Atype must be D for dichtomous or C for continuous for exposure A.")}
    if(Mtype!="C"& Mtype!="D"){stop("Error: Mtype must be D for dichtomous or C for continuous for mediator M.")}
    if(Ytype!="C"& Ytype!="D"){stop("Error: Ytype must be D for dichtomous or C for continuous for outcome Y.")}
    for(jj in 1:length(Ctype)){if(Ctype[jj]!="C" & Ctype[jj]!="D"){stop("Error: Ctype must be D for dichtomous or C for continuous measured confounders C.")}}
    for(jj in 1:length(Utype)){if(Utype[jj]!="C" & Utype[jj]!="D"){stop("Error: Utype must be D for dichtomous or C for continuous unmeasured confounders U.")}}
 
     ##################
     # check A variables
     ##################
     if((length(Atype)!=length(varA))| (length(Atype)!=length(gamma0))| (length(Atype)!=length(alphaA))|(length(Atype)!=length(betaA)) |(length(Atype)!=1 ) ){stop("Error: The length of Atype does not equal the length of varA, gamma0, alphaA, and or betaA, which does not equal one. This function does not accomodate multiple exposures A.")}
     
     ##################
     # check M variables
     ##################
     if( (length(Mtype)!=length(varM))|(length(Mtype)!=length(alpha0)) |(length(Mtype)!=length(betaM)) |(length(Mtype)!=1 ) ){stop("Error: The length of Mtype does not equal the length of varM, alpha0, and or betaM, which does not equal one. This function does not accomodate multiple mediators M.")}
    
    ##################
     # check Y variables
     ##################
     if((length(Ytype)!=length(varY))|  (length(Ytype)!=length(beta0)) |(length(Ytype)!=1)){stop("Error: The length of Ytype does not equal the length of varY, beta0, and or betaI, which does not equal one. This function does not accomodate multiple outcomes Y.")}
   if(interact==TRUE & (length(betaI)!=1)){stop("Error: length of betaI must be 1 when interact flag is TRUE.")}
   
     ##################
     # check C variables
     ##################
     if((length(Ctype)!=length(muC)) |(length(Ctype)!=length(varC)) |(length(Ctype)!=length(gammaC)) |(length(Ctype)!=length(alphaC))|(length(Ctype)!=length(betaC))){stop(paste("Error: Vectors for Ctype, muC, varC, gammaC, alphaC, and betaC must all be the same length."))}
     for(j1 in 1:length(Ctype)){
         if(Ctype[j1]=="D"){
             if(muC[j1]<0|muC[j1]==0|muC[j1]==1|muC[j1]>1){stop(paste("Error: Prob(C",j1,"=1) is ",muC[j1],". Prob(C",j1,"=1) must be between 0 and 1 for Ctype",j1,". ",sep=""))}
         }}
         
     ##################
     # check U variables
     ##################
     if((length(Utype)!=length(muU)) |(length(Utype)!=length(varU)) |(length(Utype)!=length(gammaU))|(length(Utype)!=length(alphaU))|(length(Utype)!=length(betaU))){stop(paste("Error: Vectors for Utype, muU, varU, gammaU, alphaU, and betaU must all be the same length."))}
     for(j2 in 1:length(Utype)){
         if(Utype[j2]=="D"){
             if(muU[j2]<0|muU[j2]==0|muU[j2]==1|muU[j2]>1){stop(paste("Error: Prob(U",j2,"=1) is ",muU[j2],". Prob(U",j2,"=1) must be between 0 and 1 for Utype",j2,". ",sep=""))}
         }}
     
     ##################
     #check interact flag
     ##################
     if(interact!=FALSE & interact!=TRUE){stop("Error: interact must be TRUE or FALSE. This flag is case sensitive.")}
     if(length(interact)!=1){stop("Error: length of the vector for interact must be 1.")}
 
     ##################
     #check alpha
     ##################
     if(alpha<0 | alpha>1|alpha==0|alpha==1){stop("Error: alpha must be greater than 0 and less than 1")}
     if(length(alpha)!=1){stop("Error: length of the vector for alpha must be 1")}
     
     ##################
     #check n, nSim, nBoot, seed
     ##################
     if(length(n)!=1){stop("Error: length of the vector for n must be 1")}
     if(length(nSim)!=1){stop("Error: length of the vector for nSim must be 1")}
     if(length(nBoot)!=1){stop("Error: length of the vector for nBoot must be 1")}
     if(length(seed)!=1){stop("Error: length of the vector for seed must be 1")}
     if((floor(n)!=ceiling(n))|(n<0)|(n==0)){stop(paste("Error: the sample size n must be an integer greater than 0"))}
     if((floor(nSim)!=ceiling(nSim))|(nSim<0)|(nSim==0)){stop(paste("Error: the number of simulations nSim must be an integer greater than 0"))}
     if((floor(nBoot)!=ceiling(nBoot))|(nBoot<0)|(nBoot==0)){stop(paste("Error: the number of bootstrap samples nBoot must be an integer greater than 0"))}
      if((floor(seed)!=ceiling(seed))|(seed<0)|(seed==0)){stop(paste("Error: the seed must be an integer greater than 0"))}
      
      ##################
      #atreat and acontrol
      ##################
      if(length(atreat)!=1 |length(acontrol)!=1 ){stop("Error: length of the vector for atreat and acontrol must be 1")}
      
     ###################
}#end of function


