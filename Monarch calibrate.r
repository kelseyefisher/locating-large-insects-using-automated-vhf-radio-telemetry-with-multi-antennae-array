# 072818 Monarch movement at BGE 
# estimate calibration functions, starting with raw data and key fiele

library(dplyr)
library(lattice)
library(lme4)
library(bbmle)

par(mar=c(3,3,0,0)+0.2, mgp=c(2, 0.8, 0)) 

# read data and calculate distance and angle for known locations
# jul18 <- read.csv('070718_BGE_Semicircle&RandomLocations_air&grass_4Tower4Ant_AllPower.csv', as.is=T)

jul18 <- read.csv('BGE070718.csv', as.is=T)

jul18$EstAzimuth <- jul18$EstAntAz
# copy to rename Azimuth direction

jul18<- jul18 %>% mutate(
  Distance = sqrt((T_N - L_N)^2 + (T_E - L_E)^2 + 7^2),
  logDistance = log(Distance),
  AngleD = ((EstAzimuth - CalcAzimuth + 180) %% 360) - 180,
  Angle = AngleD*2*pi/360,
  frontlobe = (AngleD > -90) & (AngleD < 90),
  backlobe = 1-frontlobe,
  anglef = frontlobe*((1-cos(Angle))^(0.95)-1),
  angleb = backlobe*((1+cos(Angle))^(0.95)-1)
  )

# distance computation assumes antennae are 7m above transmitter
#   only makes a difference for horizontal distances within ca 30m of tower

jul18locs <- jul18
# to use all data (circle, distance, and locations) for calibration

jul18.mm <- lmer(Power ~ logDistance + anglef + angleb 
  + (1|Tower:Antenna) + (0+anglef|Tower:Antenna) 
  + (0+angleb|Tower:Antenna), subset = (Distance < 150),
  data=jul18locs)

# fit linear model up to 150 m

