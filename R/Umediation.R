
Umediation<-function(
  n=100,Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",
  interact=FALSE,muC=0,varC=1,muU=0,varU=1,gamma0=0,gammaC=0,gammaU=0,varA=1,
  alpha0=0,alphaA=0,alphaC=0,alphaU=0,varM=1,
  beta0=0,betaA=0,betaM=0,betaI=0,betaC=0,betaU=0,varY=1,
  alpha=0.05,nSim=250,nBoot=400,seed=1,atreat=1,acontrol=0,
  use_multi_processing=F, num_jobs=1){
  if(num_jobs < 1){
    stop("invalid parameter: num_jobs < 1")
  }
  if(num_jobs > 1 && !use_multi_processing){
    stop("incompatible parameters: use_multi_processing is false, and num_jobs != 1")
  }
  #######################################
  # Check Input for Errors
  #######################################
  ErrorCheck(n,Atype,Mtype,Ytype,Ctype,Utype,interact,muC,varC,muU,varU,gamma0,gammaC,gammaU,varA,alpha0,alphaA,alphaC,alphaU,varM,beta0,betaA,betaM,betaI,betaC,betaU,varY,alpha,nSim,nBoot,seed,atreat,acontrol)
  
    # library(mediation) #load the mediation library
    # library(car) #load the car library

  #######################################
  # Create Results Matrix
  #######################################
  
  Results<-matrix(0,nrow=12,ncol=1)
  rownames(Results)<-c("Prop. of simulations w/ significant ACME excluding U","Prop. of simulations w/ significant ACME including U","Prop. of simulations where conclusions based on ACME match","Average ACME excluding U","Average ACME including U","Average absolute difference of ACME including U minus ACME excluding U","Prop. of simulations w/ significant ADE excluding U","Prop. of simulations w/ significant ADE including U","Prop. of simulations where conclusions based on ADE match","Average ADE excluding U","Average ADE including U","Average absolute difference of ADE including U minus ADE excluding U")
      
  #######################################
  # Create Input Data Matrix
  #######################################
  
  data_matrix_gen = generate_data_matrix(
    n=n, interact=interact,
    Atype=Atype,Mtype=Mtype,Ytype=Ytype,Ctype=Ctype,Utype=Utype,
    muC=muC,muU=muU,varC=varC,varU=varU,varM=varM,varY=varY,
    alpha0=alpha0,alphaA=alphaA,alphaC=alphaC,alphaU=alphaU,
    beta0=beta0,betaA=betaA,betaM=betaM,betaI=betaI,betaC=betaC,betaU=betaU,
    nSim=nSim,seed=seed,nBoot=nBoot
  )
  data_matrix = data_matrix_gen$data
  matC = data_matrix_gen$matC
  #######################################
  # Run the mediation and collect the results
  #######################################
  
  if(use_multi_processing){
    options(mediate.jobs = num_jobs)
    if(parallel::detectCores() == 1){
      warning("your machine may not be suitable for multiprocessing, only 1 core was detected")
    }
    if(num_jobs < 2){
      stop("There is no point in using MultiProcessing with less than 2 jobs")
    }
    if((nSim / num_jobs) < 1.0){
      warning(paste("you don't have enough Simulations in nSim:", nSim, " to fully benefit from num_jobs:", num_jobs, sep=""))
    }
    result.matrix = mediate_parallel(data_matrix)
  } else {
    if(num_jobs != 1){
      stop("num_jobs != 1 and use_multi_processing=F")
    }
    result.matrix = pbapply::pblapply(data_matrix, perform_mediation)
    dim(result.matrix) = dim(data_matrix)
  }
  rm(data_matrix)
  
  #######################################
  # Loop through the mediation result matrix
  #######################################
  
  for(si in 1:nSim){ #loop through all simulations
    
    # TODO access the med.out/ med.outU values for this simulation
    result_element = result.matrix[[si]]
    
    print(paste(is.null(result_element),names(result_element)))
    
    med.out = result_element$med.out
    med.outU = result_element$med.outU
    
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
    
  }
  
  #######################################
  # Results
  #######################################
  
  matR<-Results/nSim
  listA<-list(matR,matC,paste("Warning: correlations are only valid if at least one of the variables is normally distributed."))
  names(listA)<-c("Results","Correlations_Between_Variables","Warning")
  listA
    
} #end of function
