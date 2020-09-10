
# plots for paper

# estimated locations

#pdf(file='walk2new.pdf', height=6, width=6)

par(mfrow=c(1,1), mar=c(3,3,0,0)+0.5, mgp=c(2,0.8,0))
plot.tower(move1, label=T, center=F)
text(allwalk[,'est_E'], allwalk[,'est_N'], 1:npts, cex=0.7)
#mtext(walk, 3, 0.3, cex=0.8)

# add known locations
walk2true <- read.csv('Walk2_TrueLocations.csv', as.is=T)
points(walk2true$L_E, walk2true$L_N, pch=19, col=4)

#dev.off()


# repeat with every 15 seconds
npts <- dim(allwalk)[1]

par(mfrow=c(1,1), mar=c(3,3,0.5,0.5)+0.3, mgp=c(2,0.8,0))

plot.tower(move1, label=F, center=F)
thin <- allwalk[seq(1,npts, 3),]
thinnpts <- dim(thin)[1]
text(thin[,'est_E'], thin[,'est_N'], 1:thinnpts, cex=0.7)
#mtext(walk, 3, 0.3, cex=0.8)
#dev.off()

# add known locations
walk2true <- read.csv('Walk2_TrueLocations.csv', as.is=T)
points(walk2true$L_E, walk2true$L_N, pch=19, col=4)

# dev.off()  

# pdf(file='smoothwalk.pdf', height=5, width=7)
par(mfrow=c(1,1), mar=c(3,3,0,0)+0.3, mgp=c(2,0.8,0))

move1 <- subset(move, (L_ID==walk) )
move1$time <- hms(move1$Time)

move1$secs <- period_to_seconds(move1$time)
begin <- seconds_to_period(min(move1$secs))
end <- seconds_to_period(max(move1$secs))

bit <- subset(move1, Tower=='T1_02')
# and create "big" vector of prediction times every 5 seconds
    # in seconds

predtime <- seq(min(move1$secs), max(move1$secs), 5)
  

with(bit, plot(time, Power, type='n', pch=19) )

for (i in unique(bit$TA)) {
    bit1 <- subset(bit, TA==i)
    with(bit1, lines(time, Power, type='b', lty=3, lwd=2, pch=19) )
  
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
    print(summary(temp)$edf)
    # make predictions for all interior times
    temp.pred <- predict(temp, newdata=data.frame(secs=predtime1) )
    # and set those in large gaps to missing
    temp.pred[predgap > large] <- NA
    lines(predtime1, temp.pred, col=4, lwd=2)
    }
  }
dev.off()

# compute difference between true location and estimated locations
# match up known location times to estimated loc times

# known locations in walk2true
# estimated times in alltimes

walk2true$time <- hms(walk2true$Time)
walk2true$seconds <- period_to_seconds(walk2true$time)

matchtime <- match(walk2true$seconds, alltimes)

# double check that we have correct matching
walk2true[1:5,c('Time','L_N','L_E')]
allwalk[matchtime[1:5],]
seconds_to_period(alltimes[matchtime[1:5]])

delta <- sqrt( 
  (walk2true$L_E - allwalk[matchtime,'est_E'])^2 + 
  (walk2true$L_N - allwalk[matchtime, 'est_N'])^2 )
    
range(delta, na.rm=T)
mean(delta, na.rm=T)
median(delta, na.rm=T)

# calculate area and effective radius of 95% confint  
#matchtime <- na.omit(matchtime)
CIarea <- crit <- inout <- rep(NA, length(matchtime))

for (i in 1:length(matchtime)) {
  thisest <- allwalk[matchtime[i],]
  m <- matrix(0, nrow=2, ncol=2)
  m[1,1] <- thisest[4]
  m[1,2] <- m[2,1] <- thisest[5]
  m[2,2] <- thisest[6]
  n <- thisest[1]
  CIarea[i] <- pi*sqrt(det(m))*(2*(n-1)*qf(0.95, 2, n-2)/(n-2))
  
  # is estimate in the 95% ci?
  deltai <- cbind(walk2true$L_E - allwalk[matchtime,'est_E'],
    walk2true$L_N - allwalk[matchtime, 'est_N'] )
  crit[i] <- t(deltai[i,]) %*% solve(m) %*% deltai[i,]
  nobs <- allwalk[i,1]
  inout[i] <- (crit[i] < 2*(nobs-1)*qf(0.95, 2, nobs-2)/(nobs-2)) + 0
}

# proportion of estimates within 95% ci
sum(inout, na.rm=T)
sum(!is.na(inout))
100*mean(inout, na.rm=T)

# effective radius
effRadius <- sqrt(CIarea/pi)
median(effRadius)
range(effRadius)

# distance of each thinned estimated location from tower 1
towerE <- unique(allpred$T_E)[1]
towerN <- unique(allpred$T_N)[1]

sqrt((thin[,2] - towerE)^2 + (thin[,3] -towerN)^2)
# obs 28 is out of place: 26m from tower 1
