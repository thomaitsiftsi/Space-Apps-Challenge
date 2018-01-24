# set directory
setwd("C:/Users/Thomai/Desktop/NASA")


### load the data for pressure, humidity, solar radiation, temperature, wind direction, wind speed 

barpres <- read.table("barometric pressure.csv", header=FALSE,sep=",")

hum = read.table("humidity.csv", header=FALSE,sep=",")

solra= read.table("solar radiation.csv", header=FALSE,sep=",")

temp = read.table("temperature.csv", header=FALSE,sep=",")

winddir = read.table("wind direction in degrees.csv", header=FALSE,sep=",")

windspeed = read.table("wind speed.csv", header=FALSE,sep=",")

sunrise = read.table("sunrise.csv", header=FALSE,sep=",")

sunset = read.table("sunset.csv", header=FALSE,sep=",")


# create a data frame of the data with common unixtime, date, local time

dfmain = data.frame(barpres[,1:5], hum[,5],solra[,5], temp[,5],winddir[,5],windspeed[,5])  

dfmain = dfmain[,2:10]

dfsuntimes = data.frame(sunrise[,2:5],sunset[,5])

# create the column names for our data
colnames(dfmain) = c("unixtime", "date", "localtime", "bar. pressure", "humidity","sol. radiation","temperature","wind direction", "wind speed")

colnames(dfsuntimes) = c("unixtime", "date", "localtime","sunrise","sunset")



### find the number of observations per day


# u finds the unique dates in the column called date
u = unique(dfmain[,2])

#dim(dfmain[dfmain$date == "2016-09-29",])[1]
dfdays=c(0,0,0,0)


# calculating the daily energy outcome by finding unique days

for (dt in u)
{
  print(dt)
  x = dfmain[dfmain$date == dt,]
  avrt = 0
  tot=c(0,0,0)
  for( i in 1:(dim(x)[1]-1))
  {
    delt = 72*x[i,6] * (x[i+1,1] - x[i,1])/(3600)
    
    tot[1] = (tot[1] + delt*0.1)
    tot[2] = (tot[2] + delt*0.08)
    tot[3] = (tot[3] + delt*0.12)
    
    
  }
  
  dfdays=(rbind(dfdays, c(dt,tot[1],tot[2],tot[3])))
}

# extract the first row of zeros
dfdays = dfdays[-1,]

#plot the daily energy outcome
plot(dfdays[,2],type = "b", xlab = "Index", ylab = "Daily total energy / J", cex.lab = 1.25, col = "darkblue",bg = "lightblue", pch = 21)  




#### do EVT in the minima (which we have negated to have as maxima)

library(ismev)
library(extRemes)

# turn the energy values into numeric
dat = as.numeric(dfdays[,2])
dfweek = rep(0, floor(length(dat)/7))

#block the data into weeks and extract the weekly maxima 

for(j in 1:floor(length(dat)/7))
{
  wk = dat[(7*(j-1)):(7*j)]
  dfweek[j] = max(wk);
}


# plot maxima per week 
plot(dfweek,type = "b", xlab = "Week", ylab = "Maximum daily total energy / J", cex.lab = 1.25, col = "darkblue",bg = "lightblue", pch = 21)  



# perform the EVT analysis using package ismev


A = gev.fit(dfweek)

#diagnostic plots
gev.diag(A)




# perform EVT analysis using package extremes
fitweeks = fevd(dfweek,units = "joules")

# diagnostic plots
plot(fitweeks)

# Plots of the likelihood and its gradient for each parameter (holding the other parameters fixed
# ath their MLE estimates)
plot(fitweeks,"trace")



# parameter estimates and their CIs
ci(fitweeks,type = "parameter")



## the two-year return level corresponds to the median of the GEV
return.level(fitweeks)

# return levels and their CIs
return.level(fitweeks,do.ci = TRUE)


# return levels and their CIs for a given period of 2,3,4,5 years
ci(fitweeks,return.period = c(2,3,4,5))



#There are two ways we can test about the Gumbel hypothesis. The easiest is to use
# the following line of code which gives CI for each parameter estimate (default used below gives the 95% CI) which 
#reveals that we reject the null hypothesis of a light tailed df in favour of the upper bounded type. 




# our previous shape estimate is too close to zero so fit Gumbel and then compare models.

fitweeksgumb= fevd(dfweek,type= "Gumbel", units = "joules")


#summary of model
fitweeksgumb

#diagnostic plots
plot(fitweeksgumb)

# trace from model 
plot(fitweeksgumb,"trace")

#return levels for Gumbel model
return.level(fitweeksgumb)

#return levels for Gumbel model and their CIs
return.level(fitweeksgumb,do.ci = TRUE)


#parameter estimates for Gumbel model and their CIs
ci(fitweeksgumb,type = "parameter")

#return levels for Gumbel model for given periods of time
ci(fitweeksgumb,return.period = c(2,3,4,5))



## use likel. test to find which of the two models is favoured.Since p>a then do not reject the null (Gumbel) 
# The results of the test provide do not support for  rejecting the null hypothesis
# of the Gumbel.

lr.test(fitweeksgumb, fitweeks)


# calculate the probabilities of having energies less 

pextRemes(fitweeksgumb, q = c(-10000,0), lower.tail = FALSE)




############################## TRYING WITH EXCESSES



# turn days into numerical

days = as.numeric(dfdays[,2])


#choose threshold for our excesses
mrl.plot(days)


# confirm the above chosen threshold
threshrange.plot(days, r = c(-60000,-20000), nint = 50)



# fit a generalised Pareto distribution on the excesses
fitD = fevd(days, threshold = -35000, type = "GP", time.units = "days")


plot(fitD)

ci(fitD,type = "parameter")


ci(fitD,return.period = c(2,3,4,5))







fitexp = fevd(days,threshold = -35000, type = "Exponential", time.units = "days", verbose = TRUE)

#fit it with package ismev
test= gpd.fit(days,threshold=-35000, npy = 120, show=TRUE)
gpd.diag(test)

