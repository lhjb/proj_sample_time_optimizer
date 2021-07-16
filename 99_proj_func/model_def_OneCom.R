################################################################################
# Model definition
# one-compartment pharmacokinetic model with linear absorption 
# analytical solution, applicable for both multiple and single dosing
################################################################################

library(PopED)
## -----------------------------------------------------------------------------
# define model 
ff <- function(model_switch,xt,parameters,poped.db){
  with(as.list(parameters),{
    N = floor(xt/TAU)+1
    y=(DOSE*Favail/V)*(KA/(KA - CL/V)) * 
      (exp(-CL/V * (xt - (N - 1) * TAU)) * (1 - exp(-N * CL/V * TAU))/(1 - exp(-CL/V * TAU)) - 
         exp(-KA * (xt - (N - 1) * TAU)) * (1 - exp(-N * KA * TAU))/(1 - exp(-KA * TAU)))  
    return(list( y=y,poped.db=poped.db))
  })
}

## -----------------------------------------------------------------------------
# define the parameters
sfg <- function(x,a,bpop,b,bocc){
  # between-subject variability (BSV) for each parameter is log-normally distributed
  parameters=c( V=bpop[1]*exp(b[1]),
                KA=bpop[2]*exp(b[2]),
                CL=bpop[3]*exp(b[3]),
  # Favail is assumed not to have BSV
                Favail=bpop[4],
  # covariates (in vector a) to be optimized
                DOSE=a[1],
                TAU=a[2])
}

## -----------------------------------------------------------------------------
# define the residual unexplained variability (RUV) function
# w. both an additive and proportional component
feps <- function(model_switch,xt,parameters,epsi,poped.db){
  returnArgs <- ff(model_switch,xt,parameters,poped.db) 
  y <- returnArgs[[1]]
  poped.db <- returnArgs[[2]]
  
  y = y*(1+epsi[,1])+epsi[,2]
  
  return(list(y=y,poped.db=poped.db)) 
}


