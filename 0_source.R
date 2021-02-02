####################################################################################################
###
### File:    0_source.R
### Purpose: Load required packages and functions used throughout the SEIR model fitting
### Authors: James P. Gleeson, Thomas Brendan Murphy, Joseph D. O'Brien, and David J. P. O'Sullivan
### Date:    2/2/21
###
####################################################################################################

### packages required

packages <- c('brms', # Bayesian regression models
              'mgcv', # for GAM models
              'tidyverse', # data-wrangling + plotting
              'lubridate', # dates and time-series
              'here' # efficient file structures
              )

install.packages(setdiff(packages, rownames(installed.packages())))
lapply(packages, library, character.only = TRUE)

### functions used in SEIR model

parent_from_child <- function(C, wpc, rp, rc, delta){
  # implementation of Eq.(19) in technical note
  # C - number in the 'child' compartment
  # wpc - weight of parent-child pair
  # rp - rate associated with parent compartment
  # rc - rate associated with child compartment
  # delta - time-step used
  bigN <- length(C)
  P <- rep(0, bigN)
  for (n in 1:(bigN-1)){
    P[n] <- (1/(wpc*rp))*(((C[n+1]-C[n])/delta) + rc*C[n])
  }
  P[bigN] <- (1/(wpc*rp))*((C[bigN]-C[bigN-1])/delta + rc*C[bigN])
  P
}

child_from_parent <- function(P, wpc, rp, rc, delta){
  # implementation of Eq.(20) in technical note
  # P - number in the 'parent' compartment
  # wpc - weight of parent-child pair
  # rp - rate associated with parent compartment
  # rc - rate associated with child compartment
  # delta - time-step used
  bigN <- length(P)
  C <- rep(0, bigN)
  C0 <- 0
  C[1] <- C0 + delta*(wpc*rp*P[1] - rc*C0)
  for (n in 1:(bigN-1)){
    C[n+1] <- C[n] + delta*(wpc*rp*P[n] - rc*C[n])
  }
  C
}

## compartment change equations (Eq.(1-10))

SF_update <- function(SFn, timestep, beta, IPFn, h, IaFn, i, IiFn, It1Fn, j, It2Fn, InFn, Nv){
  # discretized Eq.(1)
  SF_new <- SFn - timestep*beta*SFn*((IpFn + h*IaFn + i*IiFn +
                                        It1Fn + j*It2Fn + InFn)/Nv)
}

EF_update <- function(EFn, timestep, beta, SFn, IPFn, h, IaFn, i, IiFn, It1Fn, j, It2Fn, InFn, Nv, L){
  # discretized Eq.(2)
  EFnew <- EFn + timestep*(beta*SFn*((IpFn + h*IaFn + i*IiFn +
                                        It1Fn + j*It2Fn + InFn)/Nv) -
                             (EFn/L))
}

IpF_update <- function(IpFn, timestep, f, L, EFn, Cv){
  # discretized Eq.(3)
  IpFnew <- IpFn + timestep*(((1-f)/L)*EFn - (IpFn/max(0.1, Cv - L)))
}

IaF_update <- function(IaFn, timestep, f, L, Efn, Dv){
  # discretized Eq.(4)
  IaFnew <- IaFn + timestep*((f/L)*Efn - (IaFn/Dv))
}

IiF_update <- function(IiFn, timestep, Cv, L, q, IpFn, Dv){
  # discretized Eq.(5)
  IiF_new <- IiFn + timestep*((1/max(0.1, Cv-L))*q*IpFn -
                                (1/max(0.1,Dv-Cv+L))*IiFn)
}

It1Fn_update <- function(It1Fn, timestep, Cv, L, tv, IpFn, T){
  # discretized Eq.(6)
  It1Fn + timestep*((1/max(0.1, Cv-L))*tv*IpFn - (It1Fn/T))
}

It2Fn_update <- function(It2Fn, timestep, It1Fn, Cv, L, Dv, T){
  # # discretized Eq.(7)
  It2Fn_new <- It2Fn + timestep*((It1Fn/T) - (It2Fn/(max(0.1, Dv - Cv + L - T))))
}

InFn_update <- function(InFn, timestep, Cv, L, tv, q, IpFn, Dv){
  # discretized Eq.(8)
  InFn_new <- InFn + timestep*(((1 - tv - q)*IpFn)/max(0.1, Cv - L) -
                                 InFn/max(0.1, Dv - Cv + L))
}

RFn_update <- function(RFn, timestep, IaFn, Dv, IiFn, Cv, L, It2Fn, T, InFn){
  # discretized Eq.(9)
  RFn_new <- RFn + timestep*((IaFn/Dv) + (IiFn/max(0.1, Dv - Cv + L))  +
                               (It2Fn/max(0.1, Dv - Cv + L - T)) +
                               (InFn/max(0.1, Dv - Cv + L)))
}

CcFn_update <- function(CcFn, timestep, It1Fn, T){
  # discretized Eq.(10)
  CcFn_new <- CcFn + timestep*(It1Fn/T);
}

### plot settings

theme_new <- function(base_size = 14, base_family = "Roboto"){
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text = element_text(size = 12, colour = "grey30"),
      legend.key=element_rect(colour=NA, fill =NA),
      axis.line = element_line(colour = 'black'),
      axis.ticks =         element_line(colour = "grey20"),
      plot.title.position = 'plot'
    )
}

