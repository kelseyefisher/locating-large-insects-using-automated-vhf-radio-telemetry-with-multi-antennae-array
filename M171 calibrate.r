# M171 calibrate.r
# code to calibrate 2019 GM field

dist <- read.csv('GM_Distance_MedianPower_Formatted.csv', as.is=T)
dist$type <- 'distance'

circ <- read.csv('GM_Semicircle_MedianPower_Formatted.csv', as.is=T)
circ$type <- 'circle'

calib <- rbind(dist, circ)

calib <- calib %>% mutate(
  TA = paste(Tower, Antenna, sep=':'),
  Distance = sqrt((T_N - L_N)^2 + (T_E - L_E)^2 + 7^2),
  logDistance = log(Distance),
  AngleD = ((EstAzimuth - CalcAzimuth + 180) %% 360) - 180,
  Angle = AngleD*2*pi/360,
  frontlobe = (AngleD > -90) & (AngleD < 90),
  backlobe = 1-frontlobe,
  anglef = frontlobe*((1-cos(Angle))^(0.95)-1),
  angleb = backlobe*((1+cos(Angle))^(0.95)-1)
  )
calibB <- filter(calib,
  Tower %in% c('T5','T6','T8','T9'))

# with(calibB, plot(Distance, Power, log='x'))
# with(subset(calibB,(type=='circle') & (Distance < 50)), plot(AngleD, Power))

GM.mm <- lmer(Power ~ logDistance + anglef + angleb 
    + (1|Tower:Antenna) + (0+anglef|Tower:Antenna) 
    + (0+angleb|Tower:Antenna), 
  subset = (Distance < 200),
  data=calib)
# fit linear model up to 200 m to all towers

