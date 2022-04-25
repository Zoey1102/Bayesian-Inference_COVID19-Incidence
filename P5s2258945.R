## Zheyi SHEN, s2258945


## ------ Overview ------
## Taking advantages of the high reliability of 
## death tolls and the fatal disease duration information,
## this project aims at making inference on new fatal/infections each day.
## More specifically, we would:
## *** Sample the hospital deaths with Covid-19 in England for 100 days since Mar.2nd, 2020
## *** Catch the posterior mean along the trajectory of "observed" new deaths n[i]
##     and also the expected new deaths m[i].
## *** Diagnostic analysis would be conducted before the final plot.
##
## Model assumption:
## The number if new infections per day changes fairly smoothly from day to day and
## the measured dist'n for fatal disease duration is correct.
## 
## Model details is given in model.jags file.


library(rjags)
library(coda)

######################################################################
######       Part I: Data reading and sampling        ################
######################################################################

## Data on daily hospital deaths with COVID-19 in England for the first 100 days from March 2nd 2020
## (Source: NHS England)
y <-c(1,2,0,2,2,0,4,4,1,9,14,20,22,27,40,46,66,63,105,103,149,159,204,263,326,353,360,437,498,576,645,647,
      700,778,743,727,813,900,792,740,779,718,699,646,686,639,610,571,524,566,486,502,451,438,386,380,345,341,
      325,313,306,268,251,259,252,266,256,215,203,196,166,183,162,180,172,167,138,160,144,153,150,122,128,116,
      133,139,122,124,117,92,83,94,109,111,83,86,83,80,73,69)

## As there's no Covid-19 deaths prior to this date, 
## we could extend the data by appending 20 days of zero deaths to the start of the data 
## New data then begins on February 11, 2020.
y_new <- c(rep(0,20),y) 
N <- length(y_new) # 120 

## Generate tmp matrix B based on the log normal density of fatal disease duration D
B <- matrix(0,N,N)  ## Initialize matrix B according to data length
## For each element B[i,j]
for (i in 1:N){  ## row i of B
    for (j in 1:N){  ## at the j-th column 
        if (j > i){  
            B[i,j] <- 0  ## 0 means no change.
        }
        else{    ## Choose the (i-j)th element in the marginal probability of log(D)
            B[i,j] <- dlnorm(i-j, meanlog = 3.235, sdlog = 0.4147, log = FALSE)
        }
    }
} ## End of B matrix construction


## Call JAGS to perform Gibbs sampling for posterior dist'n
mod <- jags.model("model.jags", data = list(y = y_new, N = N, B = B))
ns <- 10000 ## select a sample size of 10000
sam.coda <- coda.samples(mod, c("m","n"), n.iter = ns) 
## Please refer to model.jags file for model specification details


##########################################################################
######      Part II: Diagnostics        #################################
#########################################################################

## Trace plots:
par(mfrow = c(2,4), mar = c(2,2,2,2))
traceplot(sam.coda[[1]][,c("n[3]","n[40]","n[60]","n[100]", "m[3]","m[40]","m[60]","m[120]")])
## Good mixing, 
## "rapid" apparent movement about the state space could be seen. 

## ACF plots:
acfplot(sam.coda[[1]][,c("n[3]","n[40]","n[60]","n[100]", "m[3]","m[40]","m[60]","m[120]")], 
        main = "Autocorrelation Function Plots",
        aspect=2, type="l")
## The auto-correlations reduce faster around the middle part.

## For Gibbs sampling, high posterior correlation tends to slow mixing, and 
## usually the effective sample sizes are low relative to the simulation length.
## However, based on the full nodes of plot(sam.coda), the burn-in period has not ended
## (MC not yet converge) until 4000.


##  effective size of sampling calculation:
effect <- max(effectiveSize(sam.coda[[1]])); effect
## roughly 600 (it varies)



#######################################################################
##############  Part III: Summary statistics and Final Plot  ##########
#######################################################################

## Sampling for expected deaths and new infections:
death_ind <- grep("m", colnames(sam.coda[[1]])) ## index in sequential date
infection_ind <- grep("n", colnames(sam.coda[[1]]))[1:length(y)] ## using 100 here

## Credible interval for new infections:
CI <- apply(sam.coda[[1]][,infection_index], 2, quantile, prob=(c(0.025,0.975)))
upper_bound <- CI[2,]; lower_bound <- CI[1,]

## Posterior mean of new infections and expected deaths number:
expect_death <- colMeans(sam.coda[[1]][,death_ind])
new_infect <- colMeans(sam.coda[[1]][,infection_ind])

max_vertical <- max(max(new_infect), max(expect_death), 
                    max(upper_bound), max(lower_bound), 
                    max(y_new)
                    )

lockdown <- julian(as.Date("2020-3-24"), origin=as.Date("2019-12-31"))
date_vector <- seq(as.Date("2020-2-11"), by = "days", length.out = 120)

julian_day <- julian(date_vector,origin=as.Date("2019-12-31"))
infection_day <- julian_day[1:length(y)]



## The final single summary plot:
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))

plot(x=julian_day, y=y_new,
     xlab = "Day of the year",
     ylab = "Number of Individuals",
     ylim = c(0,max_vertical),
     col='grey',pch=20)

## Shade the credible interval band
polygon(x=c(rev(infection_day), infection_day), 
        y=c(rev(lower_bound), upper_bound), 
        col = "yellow", border = NA) 

lines(x=infection_day,y=lower_bound,col="black")
lines(x=infection_day,y=upper_bound,col="black")
lines(x=infection_day,y=new_infect,col="red")
lines(x=julian_day,y=expect_death,col="blue")

## Highlight Mar.24th 2020, the the first date of lock-down in the UK
abline(v=lockdown) 

text(x=lockdown + 4, y=1100,
     label="lockdown", pos=3, cex=0.6, srt=90) 


title("Daily Deaths and New Infections from Covid-19 at England (2020)", cex.main=0.8) 
legend(x="topright", legend=c("Mean of New Infections", 
                             "Mean of Expected Deaths", 
                             "95% Credible Interval (New Infections)",
                             "Observed Daily Deaths"),
       col = c("red", "blue", "yellow", "grey"), 
       lty = c(1,1,NA,NA), 
       fill = c(NA,NA, "yellow", NA),
       border = c(NA,NA,"yellow", NA), 
       pch=c(NA,NA,NA,20),
       cex = 0.55, bty = "n", x.intersp = c(2,2,1.5,2.2)
)


#######################################
### Conclusion:
###   The posterior of the new infections fit the NHS England data well.
###   As the red line is all bounded by the confidence interval.
### The tricky thing is that the peak of deaths appear after the first lockdown.
#######################################################################
