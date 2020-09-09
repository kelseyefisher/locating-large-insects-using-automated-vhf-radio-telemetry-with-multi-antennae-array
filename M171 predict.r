# M171 predict.r

# read tower information for 4 towers
M171tower <- read.csv('M171 tower.csv', as.is=T)

move1 <- read.csv('M171 hits.csv', as.is=T)
move1$time <- hms(move1$Time)

# check that every Tower/Antenna in the location data
#   matches one in M171tower
match(unique(move1$TA), unique(M171tower$TA))

# figure out beginning and end of entire run - all data
move1$secs <- period_to_seconds(move1$time)
begin <- seconds_to_period(min(move1$secs))
end <- seconds_to_period(max(move1$secs))

# data for 40 minutes, so reconstruct location every 30 secs
predtime <- seq(min(move1$secs), max(move1$secs), 30)
  
# bucket to store all predicted powers
allpred <- data.frame(Tower='A', T_N=0, T_E=0, Antenna=1, TA=1, EstAzimuth=0, 
  secs=0, Power=0, stringsAsFactors = F)

# and create "big" vector of prediction times every 30 seconds
    # in seconds

# plot power and also fit smooth curves to each time series
par(mfrow=c(2,2), mar=c(3,3,1,0)+0.2, mgp=c(2,0.8,0))
for (j in unique(move1$Tower)) {
#  j <- unique(move1$Tower)[1]
  bit <- subset(move1, Tower==j)
  with(move1, plot(time, Power, type='n', pch=19) )
  
  for (i in unique(bit$TA)) {
#  i <- unique(bit$TA)[1]
    bit1 <- subset(bit, TA==i)
    with(bit1, lines(time, Power, type='b', lty=3, lwd=2) )

    # get T_N and T_E for this tower, store in bit1
    temp <- subset(M171tower, TA==i)
    bit1$T_N <- temp$T_N
    bit1$T_E <- temp$T_E
    bit1$EstAzimuth <- temp$EstAzimuth
    
  # keep only prediction times within sampling window
    predtime1 <- subset(predtime, (predtime > min(bit1$secs)) & (predtime < max(bit1$secs)))
    npred <- length(predtime1)
    
  # find large gaps in observed data
    # define large as more than 150 seconds
    large <- 150
    start <- rev(rev(bit1$secs)[-1])   # start of each gap
    gap <- diff(bit1$secs)  # size of each gap 
    
    # determine length of gap around each prediction time
    predgap <- rep(0, npred)
    for (k in 1:npred) {
      predgap[k] <- gap[sum(start < predtime1[k])]
      }

  # fit the gam   
    if (dim(bit1)[1] > 8){
      # only fit a gam when more than 8 observed powers
    temp <- gam(Power ~ s(secs), data=bit1, gamma=2)
      # fit gam using increased smoothness (gamma > 1)
    # make predictions for all interior times
    temp.pred <- predict(temp, newdata=data.frame(secs=predtime1) )
    # and set those in large gaps to missing
    temp.pred[predgap > large] <- NA
    lines(predtime1, temp.pred, col=4)

  # concatenate this tower and antenna to output data set   
    allpred <- rbind(allpred, 
      data.frame(Tower=j, T_N=bit1$T_N[1], T_E=bit1$T_E[1], Antenna=bit1$Antenna[1], 
        TA=i, EstAzimuth=bit1$EstAzimuth[1],
        secs=predtime1, Power=temp.pred, stringsAsFactors = F))
    }
  }
  mtext(paste('Power',j), 3, 0.2, cex=0.8)
}

# and remove dummy first row
allpred <- allpred[-1,]

# how many smoothed power values for each timestamp?
table(allpred$secs)
# 5-9 values except at very end

table(allpred$secs, allpred$Tower)
# from 3 - 4 towers except at very end
#   secs with 5 obs commonly have 3 towers

# ----------------------------------

# now fit locations to smoothed power for all times with >= 5 records
M171.fit <- function(allpred, adjust) {
# fit locations using data in allpred, 
#  with adjustment to power given by adjust: 
#    0 = obs power, 22 = estimated vegetation effect

  allpredAdj <- allpred %>%
    mutate(Power = Power + adjust)
  
nobs <- table(allpredAdj$secs)
alltimes <- names(nobs)[nobs >= 5]

# bucket to keep location (E,N, and VC (as 3 element V_E, Cov, V_N))

allwalk <- matrix(NA, nrow=length(alltimes), ncol=5)
nobs <- rep(NA, length(alltimes))

for (i in 1:length(alltimes)) {
  cat(i, ' '); flush.console()
  
  bit <- subset(allpredAdj, secs == as.integer(alltimes[i]))
  bit <- na.omit(bit)  # to eliminate NA values from smoothed power

  nobs[i] <- dim(bit)[1]
  
  estloc <- optim(
    c(x=mean(bit$T_E), y=mean(bit$T_N)), 
    ss.xy, 
    hessian=T, method='BFGS', 
    ds=bit, model=GM.mm, trunc=60, s=5
    )


  allwalk[i, 1:2 ] <- estloc$par
  
  vc <- try(  solve(estloc$hessian), silent=T)

  if (class(vc) != 'try-error') {
    vc <- vc*2*estloc$value/(nobs[i]-2)   
    # adjustment to VC when using SS 
    allwalk[i, 3:5] <- c(vc[1,1], vc[1,2], vc[2,2])
    }
  }

allLocsB <- cbind(nobs, as.integer(alltimes), allwalk)

dimnames(allLocsB)[[2]] <- c('Nobs', 'Time', 'est_E', 'est_N','Var E', 'Cov','Var N')

npts <- dim(allLocsB)[1]

allLocsB
}

allLocsB0 <- M171.fit(allpred, 0)  
# estimated locations using observed power (0 adjustment)

allLocsB22 <- M171.fit(allpred, 22)
# estimated locations using veg adjustment of 22

# get Kelsey's GPSed locations
trueloc <- read.csv('TrueMonarchLocationTime_MonarchwPromisingHits_2018&2019.csv', as.is=T)
# and keep only 171
trueloc171 <- subset(trueloc, M_ID=='171')
trueloc171$secs <- period_to_seconds(hms(trueloc171$ApproximateTime))

# calculate average distance from estimate to true locations
calc.dist <- function(est, true) {
  nlocs <- dim(est)[1]
  alldist <- rep(NA, nlocs)
  
  for (i in 1:nlocs) {
    # find time matched true location
    i2 <- which.min(abs(true$secs - est[i,'Time']))
    alldist[i] <- sqrt(
      (est[i,'est_E'] - true$Easting[i2])^2 + 
      (est[i,'est_N'] - true$Northing[i2])^2)
    }
  alldist
  }


# distances using observed power
alldist0 <- calc.dist(allLocsB0, trueloc171)

# distances with 22 added
alldist22 <- calc.dist(allLocsB22, trueloc171)

