####################################################################################################
###
### File:    2_Bayesian_GAM_plots.R
### Purpose: Using the GAM SEIR model fits from Bayesian_GAM_fits.R to plot curves for
###          the dynamics of each compartment and the number of cases
### Authors: James P. Gleeson, Thomas Brendan Murphy, Joseph D. O'Brien, and David J. P. O'Sullivan
### Date:    2/2/21
###
####################################################################################################

source('0_source.R')

dir.create(file.path(here('Plots',file_date)))

apply(dailyNewCases_samples,1,quantile,prob=c(0.025,0.25,0.5,0.75,0.975)) %>%
  t() %>%
  as_tibble() %>%
  mutate(time = row_number(),
         average = dailyNewCases_samples %>% t() %>% colMeans()) %>%
  rename('q025' = `2.5%`,
           'q25' = `25%`,
           'q50' = `50%`,
           'q75' = `75%`,
           'q975' = `97.5%`) %>%
  ggplot(aes(x = time, y = q50)) +
  geom_point(data = case_data,
            aes(x = days_in, y = daycount),
            size = 0.6) +
  geom_line(color = '#ED1C24',
            size = 0.5) +
  geom_ribbon(aes(ymin = q25, ymax = q75),
              fill = '#ED1C24',
              alpha = 0.4) +
  geom_ribbon(aes(ymin = q025, ymax = q975),
              fill = '#ED1C24',
              alpha = 0.2) +
  labs(x = 'Time (Days)',
       y = NULL,
       title = 'Daily new cases as a function of days from Feb 28th') +
  theme_new()

ggsave(here('Plots',file_date, paste('daily_cases_', file_date, '.pdf', sep = '')),
       width = 6, height = 4,
       units = 'in', device = cairo_pdf)
       
dir_in <- here('Output_data/sample_curves', file_date, '/')
files <- list.files(dir_in)
plot_names <- c('Beta', 'CcF', 'Exposed',
                'Force of Infection', 'Infected',
                'Daily New Cases', 'Recovered',
                'R(t)', 'Susceptible')

for(ii in 1:length(files)){
  file <- files[ii]
  variable <- str_remove(file,'_samples.csv')
  data <- read.csv(paste0(dir_in,'/',file))
  
  apply(data,1,quantile,prob=c(0.025,0.25,0.5,0.75,0.975)) %>%
    t() %>%
    as_tibble() %>%
    mutate(time = row_number(),
           average = data %>% t() %>% colMeans()) %>%
    rename('q025' = `2.5%`,
           'q25' = `25%`,
           'q50' = `50%`,
           'q75' = `75%`,
           'q975' = `97.5%`) %>%
    ggplot(aes(x = time, y = average)) +
    geom_line(color = '#ED1C24',
              size = 0.75) +
    geom_ribbon(aes(ymin = q25, ymax = q75),
                fill = '#ED1C24',
                alpha = 0.4) +
    geom_ribbon(aes(ymin = q025, ymax = q975),
                fill = '#ED1C24',
                alpha = 0.2) +
    labs(x = 'Time (Days)',
         y = NULL,
         title = paste(plot_names[ii], 'as a function of days from Feb 28th', sep = ' ')) +
    theme_new()
  
  ggsave(here('Plots',file_date,paste(variable,'.pdf', sep ='')), width = 6, height = 4,
         units = 'in', device = cairo_pdf)
}

