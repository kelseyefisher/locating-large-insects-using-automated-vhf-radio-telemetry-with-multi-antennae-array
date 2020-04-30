# Monarch predict.r

# Estimate stationary monarch (070718 data) using median power
#   and tweak power in 3 ways
#   add 2, add 5, add data for tower 3 (SE corner)

Mtrue <- data.frame(L_E=612173.333, L_N=4667941.353)
# true location, from Kelsey

Mmove1Median <- read.csv('070718_Power_MonarchTime_M100_Formatted_MedianPower.csv',
  as.is=T)
Mmove1Median$EstAzimuth <- Mmove1Median$EstAntAzimuth
# median power for M100, stationary

bit1 <- bit2 <- bit3 <- Mmove1Median

# add 9.6 to power and call it location 102
bit1$Power <- bit1$Power + 9.6
bit1$L_ID <- 102

# add 5 to power and call it location 105
bit2$Power <- bit2$Power + 5
bit2$L_ID <- 105

# figure out missing tower information
# jul18 has calibration information - created by calibrate.r

# check that tower locations T1, T2 and T4 match in calib and monarch data
temp <- split(jul18, jul18$Tower)
temp[[3]][1,]
# they do

# copy T1 data and rename as T3
T3 <- bit3[bit3$Tower=='T1',]
T3$Tower='T3'
T3$Antenna=2
T3$T_Northing=temp[[3]]$T_N[1]
T3$T_Easting=temp[[3]]$T_E[1]
T3$EstAzimuth=70 + 180
bit3$L_ID = 200
T3$L_ID = 200

Mmove1Median <- rbind(
  Mmove1Median, bit1, bit2, bit3, T3)

alllocs <- unique(Mmove1Median$L_ID)
nlocs <- length(alllocs)
estlocMedian <- NULL

for (i in 1:nlocs) {
  bit <- subset(Mmove1Median, L_ID==alllocs[i])
  
  cat(i, ' '); flush.console()
  
  temp <-  optim(
    c(x=mean(bit$T_E), y=mean(bit$T_N)), 
    ss.xy, 
    hessian=T, method='BFGS', 
    ds=bit, model=jul18.mm, trunc=54, s=5
    )
  n <- dim(bit)[1]
  temp$vc <- solve(temp$hessian)*2*temp$value/(n-2)

  estlocMedian <- rbind(
    estlocMedian, 
    c(temp$par, 
      temp$vc[c(1,2,4)],
      dim(bit)[1])
    )
}

