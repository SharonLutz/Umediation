\name{Umediation}
\alias{Umediation}
\title{Umediation}

\description{
    The Umediation function examines the role of unmeasured confounding on the estimates for the average causal mediated effect (ACME) and average direct effect (ADE) in mediation analysis. User input specifies the relationship between the exposure A, the mediator M, the outcome Y, the measured confounder C, and the unmeasured confounder U. The function runs mediation analysis including and excluding the unmeasured confounder U in the model. Umediation allows the user to examine how the results of the mediation analysis would change if the unmeasured confounder U was included or not included in the model. The function allows for continuous or dichotomous exposure A, mediator M, outcome Y, measured confounder C, and unmeasured confounder U. Umediation allows for multiple measured confounders C and unmeasured confounders U. In addition, Umediation allows for an interaction between the exposure A and the mediator M on the outcome Y.
}
\usage{
    Umediation<-function (n=100,Atype="D",Mtype="C",Ytype="C",Ctype="C",Utype="C",interact=FALSE,muC=0,varC=1,muU=0,varU=1, 
                          gamma0=0,gammaC=0,gammaU=0,varA=1,alpha0=0,alphaA=0,alphaC=0,alphaU=0,varM=1,beta0=0,betaA=0,betaM=0,
                          betaI=0,betaC=0,betaU=0,varY=1,alpha=0.05,nSim=250,nBoot=400)
}
\arguments{
    \item{n}{
        is the sample size of the population that is being simulated.
    }
    \item{Atype}{
        is either "C" for continuous, normally distributed exposure A or "D" for dichotomous, binary exposure A (i.e. Atype=”D”).
    }
    \item{Mtype}{
        is either "C" continuous, normally distributed mediator M or "D" for dichotomous, binary mediator M (i.e. Mtype=”C”).
    }
\item{Ytype}{
        is either "C" for continuous, normally distributed outcome Y or "D" for dichotomous, binary outcome Y (i.e. Ytype=”D”).
    }
\item{Ctype}{
        is either "C" for continuous, normally distributed measured confounder C or "D" for dichotomous, binary measured confounder C. Ctype can be a single value (i.e. Ctype=”C”) or a vector for multiple measured confounders (i.e. Ctype=c(”C”,”C”,”D”) ).
    }
\item{Utype}{
is either "C" for continuous, normally distributed unmeasured confounder U or "D" for dichotomous, binary unmeasured confounder U. Utype can be a single value (i.e. Utype=”C”) or a vector for multiple measured confounders (i.e. Utype=c(”C”,”C”,”D”) ).
    }
    \item{interact}{
        Using the flag interact=TRUE allows for an interaction between the exposure A and the mediator M on the outcome Y (i.e. E[Y]=beta0+betaA*A+betaM*M+betaC*C+betaU*U+betaI*A*M). By default, interact=FALSE (i.e. E[Y]=beta0+betaA*A+betaM*M+betaC*C+betaU*U).
    }
\item{muC}{
is the mean vector for the measured confounder C. For continuous measured confounder (i.e. Ctype=”C”), muC is the mean of C. For dichotomous measured confounder C (i.e. Ctype=”D”), muC is the probability C=1. 
}
\item{varC}{
is the variance of the measured confounder C when Ctype=”C”. For multiple measured confounders, the length of varC must match muC and Ctype (i.e. Ctype=c(”C”,”C”,”D”) and muC=c(-0.1,0.2,0.3) and varC=c(1,1,1))

}
\item{muU}{
is the mean vector for the unmeasured confounder U. For continuous unmeasured confounder (i.e. Utype=”C”), muU is the mean of U. For dichotomous unmeasured confounder U (i.e. Utype=”D”), muU is the probability U=1. 
}
\item{varU}{
is the variance of the unmeasured confounder U when Utype=”C”. For multiple unmeasured confounders, the length of varU must match muU and Utype (i.e. Utype=c(”C”,”C”,”D”) and muU=c(-0.1,0.2,0.3) and varU=c(1,1,1))
}
    \item{gamma0}{
       specifies the intercept for the exposure A (i.e. logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C).
    }
    \item{gammaC}{
        specifies the relationship between the measured confounder C and the exposure A (i.e. logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C).
    }
    \item{gammaU}{
        specifies the relationship between the measured confounder U and the exposure A (i.e. logit(P(A=1)) or E[A]=gamma0+gammaU*U+gammaC*C).
    }
    \item{varA}{
        is the variance of the exposure A when Atype=”C”. Default is varA=1.
    }
    \item{alpha0}{
        specifies the intercept for the mediator M (i.e. Logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U).
    }
    \item{alphaA}{
        specifies the relationship between the exposure A and the mediator M (i.e  logit(P(M=1)) or E[M]=alpha0+alphaA*A+alphaC*C+alphaU*U).
    }
    \item{alphaC}{
       specifies the relationship between the measured confounder C and the mediator M (i.e. E[M] or logit(P(M=1))=alpha0+alphaA*A+alphaC*C+alphaU*U).
    }
    \item{alphaU}{
        specifies the relationship between the unmeasured confounder U and the mediator M (i.e. E[M] or logit(P(M=1))=alpha0+alphaA*A+alphaC*C+alphaU*U).
    }
    \item{varM}{
        is the variance of the mediator M when Mtype=”C”. Default is varM=1.
    }
    \item{beta0}{
        specifies the intercept for the outcome Y (i.e. logit(P(Y=1)) or E[Y]=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
    }
    \item{betaA}{
        specifies the relationship between the exposure A and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
    }
    \item{betaM}{
        specifies the relationship between the mediator M and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
    }
    \item{betaI}{
        specifies the interaction between the mediator M and the exposure A on the outcome Y (i.e. E[Y] or logit(P(Y=1))= beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U) when the flag interact=TRUE.
    }
    \item{betaC}{
       specifies the relationship between the measured confounder C and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
    }
    \item{betaU}{
        specifies the relationship between the unmeasured confounder U and the outcome Y (i.e. E[Y] or logit(P(Y=1))=beta0+betaA*A+betaM*M+betaI*A*M+betaC*C+betaU*U).
    }
    \item{varY}{
        is the variance of the outcome Y when Ytype=”C”. Default is varY=1.
    }
\item{alpha}{
        is the significance level. Default value is alpha=0.05.
    }
    \item{nSim}{
        is the number of simulations run for the function. The more simulations run, the more accurate the results, but this will make the function slower.
    }
    \item{nBoot}{
       is the number of Monte Carlo draws for nonparametric bootstrap or quasi-Bayesian approximation for the mediate function.
    }
    
    \item{seed}{
       sets the seed used for the random generator.
    }
    \item{atreat}{
       sets the treatment group for the exposure A.
    }
    \item{acontrol}{
       sets the control group for the exposure A.
    }
    
    
}

\value{
    The function outputs (1) the proportion of simulations where the average causal mediation effect (ACME) is significant when the model does NOT include U, (2) the proportion of simulations where the ACME is significant when the model includes U, and (3) the proportion of simulations where conclusions based on the ACME match (i.e. the ACME is significant when U is excluded from the model and included in the model or the ACME is not significant when U is excluded from the model and included in the model). The function also outputs (1) the average estimate of the average ACME when U is NOT included in the model, (2) the average ACME when U is included in the model, and (3) the average absolute difference for the ACME when U is included in the model and the ACME when U is excluded from the model. This is given for both the ACME and the average direct effect (ADE). The correlation between variables is also given to show how the change in betas, alphas, and gammas effect the relationship between these variables. Note: correlation is valid if at least on of the variables is normally distributed.
}

\author{
    Sharon Lutz, Annie Thwing
}

\section{Warning }{
    library(mediation) and library(car) are needed to run this function. 
}

\examples{
    testM<- Umediation(n=1000,Atype="D",Mtype="C",Ytype="C",Ctype=c("C","D","D"),Utype=c("C","D","D","C"),interact=TRUE,
                       muC=c(0.1,0.3,0.2),varC=c(1,1,1),muU=c(.1,0.3,0.2,.1),varU=c(1,1,1,1),gamma0=0,gammaC=c(1,0.3,0.2),
                       gammaU=c(1,0.3,0.2,0.4),varA=1,alpha0=0,alphaA=1,alphaC=c(0.3,0.2,0.2),alphaU=c(0.3,0.2,0.3,0.2),
                       varM=1,beta0=0,betaA=-1,betaM=1,betaI=1,betaC=c(0.3,0.2,0.1),betaU=c(0.3,0.2,-1.3,0.2),varY=1,alpha=0.05,
                       nSim=100,nBoot=400)
    testM
}

\keyword{ function }
\keyword{ mediation }
\keyword{ unmeasured } 
\keyword{confounding }
