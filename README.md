# Bayesian-Inference_COVID19-Incidence
Bayesian Inference on Covid-19 incidence in England

# Model Idea
Data on deaths each day from Covid-19 provide some of the most reliable data on the state of the epidemic, because the data constitute a fairly complete record of fatal cases. There is also reliable information available on the distribution of time from infection to death in fatal cases. This opens up the possibility of inferring the number of new (fatal) infections each day, from the data on fatalities each day. This can not be done in an entirely model free way, but it is sensible to use models that attempt to make the minimum of assumptions.

One general approach is to base a model only on the assumptions that the number of new infections per day changes fairly smoothly from day to day, and that the measured distribution of time from infection to death is correct.  Because of the long disease durations, using the density as a model for the probability of a disease duration rounded to the nearest day is reasonable. 

Note that while the same analysis can be conducted including deaths outside hospital, the interpretation is then more problematic: deaths outside hospital are usually of very frail patients in care homes, where the infection process does not reflect the general community transmission of most interest, and where the hospital record derived fatal disease duration distributions may not be applicable.

# Data and Model specification 
See the task description in the pdf file PS5.pdf

# Files
1. model.jags: the JAGS model specification file
2. P5s2258945.R: the R file containing the code to run my JAGS model and produce the specified results plot
