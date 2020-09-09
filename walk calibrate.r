# estimate calibration functions, starting with raw data

# read data and calculate distance and angle for known locations

apr19 <- read.csv('050319_Spring19Callibration_AllDays_MedianPower_Formatted.csv', as.is=T)
apr19key <- read.csv('051319_key.csv', as.is=T)

# The key file specifies the type of location
#   circle: different angles, 25m from tower
#   distance: one angle, different distances from tower
#   location: random locations

apr19b <- left_join(apr19, apr19key, by='L_ID')
apr19b <- apr19b %>% mutate(
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

apr19locs <- subset(apr19b, type %in% c('circle','distance'))
# use circle and distance values for calibration

apr19.mm <- lmer(Power ~ logDistance + anglef + angleb 
  + (1|Tower:Antenna) + (0+anglef|Tower:Antenna) 
  + (0+angleb|Tower:Antenna), subset = (Distance < 150),
  data=apr19locs)
# fit linear model up to 150 m

save(apr19.mm, file='apr19Model.sav')


