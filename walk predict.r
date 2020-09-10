move <- read.csv('050619_Spring19Callibration_MovementWalks_AllPower_Formatted.csv', as.is=T)
table(move$L_ID)

move$TA <- paste(move$Tower, move$Antenna, sep=':')

walk <- 'M2_02_01'
# specify which walk to fit

move1 <- subset(move, L_ID==walk)
move1$time <- hms(move1$Time)

# figure out beginning and end of entire run - all data
move1$secs <- period_to_seconds(move1$time)
begin <- seconds_to_period(min(move1$secs))
end <- seconds_to_period(max(move1$secs))

# bucket to store all predicted powers
allpred <- data.frame(Tower='A', T_N=0, T_E=0, Antenna=1, TA=1, EstAzimuth=0, 
  secs=0, Power=0, stringsAsFactors = F)

# and create "big" vector of prediction times every 5 seconds
    # in seconds

predtime <- seq(min(move1$secs), max(move1$secs), 5)
  
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
  
  # keep only prediction times within sampling window
    predtime1 <- subset(predtime, (predtime > min(bit1$secs)) & (predtime < max(bit1$secs)))
    npred <- length(predtime1)
    
  # find large gaps in observed data
    # define large as more than 60 seconds
    large <- 60
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

# now fit locations to smoothed power for all times with > 5 records

nobs <- table(allpred$secs)
alltimes <- names(nobs)[nobs > 5]

# bucket to keep location (E,N, and VC (as 3 element V_E, Cov, V_N))

allwalk <- allwalkR <- matrix(NA, nrow=length(alltimes), ncol=5)
nobs <- rep(NA, length(alltimes))

for (i in 1:length(alltimes)) {
  cat(i, ' '); flush.console()
  
  bit <- subset(allpred, secs == as.integer(alltimes[i]))
  bit <- na.omit(bit)  # to eliminate NA values from smoothed power

  nobs[i] <- dim(bit)[1]
  
  estloc <- optim(
    c(x=mean(bit$T_E), y=mean(bit$T_N)), 
    ss.xy, 
    hessian=T, method='BFGS', 
    ds=bit, model=apr19.mm, trunc=54, s=5
    )


  allwalk[i, 1:2 ] <- estloc$par
  
  vc <- try(  solve(estloc$hessian), silent=T)

  if (class(vc) != 'try-error') {
    vc <- vc*2*estloc$value/(nobs[i]-2)   
    # adjustment to VC when using SS 
    allwalk[i, 3:5] <- c(vc[1,1], vc[1,2], vc[2,2])
    }
  }

allwalk <- cbind(nobs, allwalk)

dimnames(allwalk)[[2]] <- c('Nobs', 'est_E', 'est_N','Var E', 'Cov','Var N')

npts <- dim(allwalkR)[1]

