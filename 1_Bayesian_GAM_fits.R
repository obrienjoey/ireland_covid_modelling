####################################################################################################
###
### File:    1_Bayesian_GAM_fits.R
### Purpose: Using the daily case data, fit the GAM via a Bayesian framework as outlined in
###          'A population-level SEIR model for COVID-19 scenarios (updated)
###          returning the dynamics of each compartment over the time window of interest.
###          Furthermore, can forecast scenarios via an reproductive number of choice from
###          today's date
### Authors: James P. Gleeson, Thomas Brendan Murphy, Joseph D. O'Brien, and David J. P. O'Sullivan
### Date:    2/2/21
###
####################################################################################################

source('0_source.R')

### parameters

steps_in_day <- 10 # no. timesteps in one day (INTEGER)
timestep <- 1/steps_in_day # time step at a daily level
no_knots <- 10 # number of knots for the Bayesian GAM

### actual case data

file_in <- 'case_data_IRL_20201111.csv' # filename must have this format
file_date <- parse_number(file_in)
case_data <- read_csv(here('Input_data', file_in))

### set up the data

case_data <- case_data %>%
  mutate(days_in = seq(0,nrow(.)-1))

### Fit the Bayesian model as outlined in Sec.2

fits_br <- brm(data = case_data,
               family = negbinomial(),
               daycount~s(days_in,bs="tp",k=no_knots),
               iter=4000,
               control=list(adapt_delta=0.99,max_treedepth=20))

### save the model fit
filename <- here('Fits',paste0("brmsfit","_",file_date,".Rdata"))
save(fits_br,file=filename)

### Generate the case curves from the posterior distribution
daygrid <- seq(0,nrow(case_data)-1,by=timestep)
predmat <- posterior_epred(fits_br,newdata=list(days_in=daygrid))

no_samples <- 1000 # number of draws from posterior # (JOB: in the note this is 10000?)
ind <- sample(1:nrow(predmat),size=no_samples)
paths <- t(predmat)[,ind] # samples drawn
nt <- length(daygrid) # total number of time-steps
endT <- floor(timestep*nt) # final time

### Scenario parameters
scenarioR0 <- 0.5 # assume Rt takes this value from today
endTscenario <- endT + 12*7 # run scenario for 12 weeks beyond today

### Now create matrices to hold the daily counts in each case
dailyCcF_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyS_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyE_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyR_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyBeta_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyI_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyFI_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyRt_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)
dailyNewCases_samples <- matrix(0, nrow = endTscenario, ncol = no_samples)

### now run the simulation across each sample from the posterior
for(sample in 1:no_samples){ # for each sample curve
  gv <- paths[,sample]
  
  ### Draw the parameters from a uniform distribution as per Sec.1.1
  L <- runif(1,3.9,5.9) # average latent period
  Cv <- runif(1,max(L,4.8), 6.8) # average incubation period
  Dv <- runif(1,max(Cv-L,5),9) # average infectious period
  h <- runif(1,0.01,0.5) # multip. factor for reduction of eff. trans. from the
  # Asymp. I compartment, relative to Sympt. I
  i <- runif(1,0,0.1) # multip. factor for reduction of eff. trans. from the
  # Immediate isolation compartment, relative to Sympt. I
  j <- runif(1,0,0.1) # multip. factor for reduction of eff. trans. from the
  # Post-test isolation compartment, relative to Sympt. I
  f <- runif(1,0.18,0.82) # fraction of infected who are Asymptomatic
  tv <- runif(1,0.5,1) # frac of sympt. cases that are tested
  q <- runif(1,0,1-tv) # frac of sympt. cases that self-quar. from app. of symptoms
  # until recovery
  T <- runif(1,1,max(5,Dv-Cv+L)) # average wait for test results
  Nv <- 4.9e6 # population
  
  R0_over_beta <- (f-1)*((i-1)*(Cv-L)*q + (j-1)*(Cv-L+T)*tv) +
    Dv*(1 + (i-1)*q + (j-1)*tv + f*(-1+h+q-(i*q)+tv-(j*tv))) # Eq.(22)
  
  ### Invert SEIR equations to obtain time-varying beta(t) as per Sec.2.3
  It1 <- T*gv # from eq.10 of the Technical Note
  IpPC <- parent_from_child(It1, tv, 1/max(0.1,(Cv-L)), 1/T, timestep)
  EPC <- parent_from_child(IpPC, 1-f, 1/L, 1/max(0.1,(Cv-L)), timestep)
  wPC <- parent_from_child(EPC, 1, 1, 1/L, timestep)
  SPC <- Nv - timestep*cumsum(wPC) # integrate -wPC to get SPC, with initial value of Nv
  
  It2PC <- child_from_parent(It1, 1, 1/T, 1/max(0.1,(Dv-Cv+L-T)), timestep)
  IiPC <- child_from_parent(IpPC, q, 1/max(0.1,(Cv-L)), 1/max(0.1,(Dv-Cv+L)), timestep)
  InPC <- child_from_parent(IpPC, 1-q-tv, 1/max(0.1,(Cv-L)), 1/max(0.1,(Dv-Cv+L)), timestep)
  IaPC <- child_from_parent(EPC, f, 1/L, 1/Dv, timestep)
  
  betavalsPC <- (-Nv)*(-wPC/SPC)*(1/(IpPC + h*IaPC + i*IiPC + It1 + j*It2PC + InPC)) # Eq.(16)
  R0valsPC <- R0_over_beta*betavalsPC
  
  # set up vectors to hold daily values
  dailyCcF <- c(rep(0,endTscenario))
  dailyS <- c(rep(0,endTscenario))
  dailyE <- c(rep(0,endTscenario))
  dailyR <- c(rep(0,endTscenario))
  dailyBeta <- c(rep(0,endTscenario))
  dailyI <- c(rep(0,endTscenario))
  
  # initial values for forward solve
  SFn <- SPC[1]
  EFn <- EPC[1]
  IpFn <- IpPC[1]
  IaFn <- IaPC[1]
  IiFn <- IiPC[1]
  It1Fn <- It1[1]
  It2Fn <- It2PC[1]
  InFn <- InPC[1]
  RFn <- 0
  CcFn <- 0
  
  n <- 0
  # solve forward from time 0 to today
  for(day in 1:(endT)){# loop over days
    for(day_time_step in 1:steps_in_day){# loop over timesteps in days
      n <- n + 1
      # use finite-difference equations
      SFn_new <- SF_update(SFn, timestep, betavalsPC[n], IPFn, h, IaFn, i, IiFn, It1Fn, j, It2Fn, InFn, Nv)
      EFn_new <- EF_update(EFn, timestep, betavalsPC[n], SFn, IPFn, h, IaFn, i, IiFn, It1Fn, j, It2Fn, InFn, Nv, L)
      IpFn_new <- IpF_update(IpFn, timestep, f, L, EFn, Cv)
      IaFn_new <- IaF_update(IaFn, timestep, f, L, EFn, Dv)
      IiFn_new <- IiF_update(IiFn, timestep, Cv, L, q, IpFn, Dv)
      It1Fn_new <- It1Fn_update(It1Fn, timestep, Cv, L, tv, IpFn, T)
      It2Fn_new <- It2Fn_update(It2Fn, timestep, It1Fn, Cv, L, Dv, T)
      InFn_new <- InFn_update(InFn, timestep, Cv, L, tv, q, IpFn, Dv)
      RFn_new <- RFn_update(RFn, timestep, IaFn, Dv, IiFn, Cv, L, It2Fn, T, InFn)
      CcFn_new <- CcFn_update(CcFn, timestep, It1Fn, T)
      
      ### update the values
      SFn <- SFn_new
      EFn <- EFn_new
      IpFn <- IpFn_new
      IaFn <- IaFn_new
      IiFn <- IiFn_new
      It1Fn <- It1Fn_new
      It2Fn <- It2Fn_new
      InFn <- InFn_new
      RFn <- RFn_new
      CcFn <- CcFn_new
    }
    dailyCcF[day] <- CcFn
    dailyS[day] <- SFn
    dailyE[day] <- EFn
    dailyR[day] <- RFn
    dailyBeta[day] <- betavalsPC[n]
    dailyI[day] <- IpFn + IaFn + IiFn + It1Fn + It2Fn + InFn
  }
  
  ### NB - Can change timestep for scenario solve here
  # uncomment following two lines and edit
  #timestep = 0.0001
  #steps_in_day=floor(1/timestep)
  # Solve forward from today according to scenario specified
  # Solver is identical to above, except beta is specified by the
  # desired value of reproductive number, see scenarioBeta usage below
  scenarioBeta <- scenarioR0/R0_over_beta
  
  for(day in (endT+1):(endTscenario)){# loop over days
    for(day_time_step in 1:steps_in_day){# loop over timesteps in days
      n <- n + 1
      # use finite-difference equations
      SFn_new <- SF_update(SFn, timestep, scenarioBeta, IPFn, h, IaFn, i, IiFn, It1Fn, j, It2Fn, InFn, Nv)
      EFn_new <- EF_update(EFn, timestep, scenarioBeta, SFn, IPFn, h, IaFn, i, IiFn, It1Fn, j, It2Fn, InFn, Nv, L)
      IpFn_new <- IpF_update(IpFn, timestep, f, L, EFn, Cv)
      IaFn_new <- IaF_update(IaFn, timestep, f, L, EFn, Dv)
      IiFn_new <- IiF_update(IiFn, timestep, Cv, L, q, IpFn, Dv)
      It1Fn_new <- It1Fn_update(It1Fn, timestep, Cv, L, tv, IpFn, T)
      It2Fn_new <- It2Fn_update(It2Fn, timestep, It1Fn, Cv, L, Dv, T)
      InFn_new <- InFn_update(InFn, timestep, Cv, L, tv, q, IpFn, Dv)
      RFn_new <- RFn_update(RFn, timestep, IaFn, Dv, IiFn, Cv, L, It2Fn, T, InFn)
      CcFn_new <- CcFn_update(CcFn, timestep, It1Fn, T)
      
      ### update the values
      SFn <- SFn_new
      EFn <- EFn_new
      IpFn <- IpFn_new
      IaFn <- IaFn_new
      IiFn <- IiFn_new
      It1Fn <- It1Fn_new
      It2Fn <- It2Fn_new
      InFn <- InFn_new
      RFn <- RFn_new
      CcFn <- CcFn_new
    }
    # store the end-of-day values for output
    dailyCcF[day] <- CcFn
    dailyS[day] <- SFn
    dailyE[day] <- EFn
    dailyR[day] <- RFn
    dailyBeta[day] <- scenarioBeta
    dailyI[day] <- IpFn + IaFn + IiFn + It1Fn + It2Fn + InFn
  }
  
  dailyFI <- dailyBeta*dailyI/Nv # calculate force of infection
  dailyRt <- R0_over_beta*dailyBeta # translate daily beta values to daily Rt values
  dailyNewCases <- diff(c(0,dailyCcF)) # difference the accumulated cases to find the new  cases each day
  
  ### store the results
  
  dailyCcF_samples[,sample] <- dailyCcF
  dailyS_samples[,sample] <- dailyS
  dailyE_samples[,sample] <- dailyE
  dailyR_samples[,sample] <- dailyR
  dailyBeta_samples[,sample] <- dailyBeta
  dailyI_samples[,sample] <- dailyI
  dailyFI_samples[,sample] <- dailyFI
  dailyRt_samples[,sample] <- dailyRt
  dailyNewCases_samples[,sample] <- dailyNewCases
}

### save the curves
dir.create(file.path(here('Output_data/sample_curves', file_date)), showWarnings = FALSE)
data_dir <- here('Output_data/sample_curves', file_date)

curve_list <- list(dailyCcF_samples,dailyS_samples,dailyE_samples,
                   dailyR_samples,dailyBeta_samples,dailyI_samples,
                   dailyFI_samples,dailyRt_samples,dailyNewCases_samples)
fnames <- c('dailyCcF_samples','dailyS_samples','dailyE_samples',
            'dailyR_samples','dailyBeta_samples','dailyI_samples',
            'dailyFI_samples','dailyRt_samples','dailyNewCases_samples')

for(file in 1:length(fnames)){
  write.csv(curve_list[[file]],
            file = paste0(data_dir,'/',fnames[file],'.csv'),
            row.names=FALSE)
}
