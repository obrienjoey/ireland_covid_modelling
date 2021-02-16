# ireland_covid_modelling

2/2/21 - Code associated with [*A population-level SEIR model for COVID-19 scenarios (updated)*](https://assets.gov.ie/74595/e10ea7bb423c4110839e2f30c3817dc7.pdf) by James P. Gleeson, Thomas Brendan Murphy, Joseph D. O'Brien, and David J. P. O'Sullivan.

Using daily Covid-19 case data, the code fits a generalized additive model (GAM) via a Bayesian framework before running simulations of the population-level susceptible-exposed-infected-removed (SEIR) model that is used by the Irish Epidemiological Modelling Advisory Group (IEMAG) reporting to the National Public Health Emergency Team (NPHET).

There are three scripts which must be run in order

1. 0_source.R - loads required libraries and defines the functions used throughout the analysis.
2. 1_Bayesian_GAM_fits.R - fits the GAM to a provided dataframe of case date via a Bayesian framework before simulating 1000 realizations of the SEIR process using draws from    the posterior distribution. These simulations are run for a specified number of days in the future to generate the future sizes of the compartments in the SEIR model, assuming that the reproduction number remains constant at the value given by parameter scenarioR0.
3. 2_Bayesian_GAM_plots.R - returns the plots of the resulting time series.

The code runs a scenario from the date 11/11/2020, until a specified time in the future. The input data supplied is a daily time-series or the number of national cases up until that point.  
