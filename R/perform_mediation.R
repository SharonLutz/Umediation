perform_mediation <- function(in_data){
  Mtype = in_data$Mtype
  Ytype = in_data$Ytype
  interact = in_data$interact
  
  M = in_data$M
  Y = in_data$Y
  A = in_data$A
  CC = in_data$CC
  U = in_data$U
  nBoot = in_data$nBoot
  
  old_rand_state <- NULL
  
  genv = globalenv()
  
  if(exists(".Random.seed", envir = genv)){
    old_rand_state <- genv[[".Random.seed"]]
    # print("storing old rand state")
    # print(old_rand_state)
  }
  
  genv[[".Random.seed"]] <- in_data$RandStateParam
  
  # print(genv[[".Random.seed"]])
  
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
  
  med.out <- mediation::mediate(med.fit, out.fit, treat = "A", mediator = "M", sim=nBoot)
  
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
  
  med.outU <- mediation::mediate(med.fitU, out.fitU, treat = "A", mediator = "M", sim=nBoot)
 
  
  if(!is.null(old_rand_state)){
    # print("restoring old rand state")
    # print(old_rand_state)
    
    genv[[".Random.seed"]] <- old_rand_state
    # print(genv[[".Random.seed"]])
  }
   
  return(list(med.out=med.out, med.outU=med.outU))
  
}