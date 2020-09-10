# new code for BGE

bge <- read.csv('BGE070718.csv', as.is=T)

bge$EstAzimuth <- bge$EstAntAz
# copy to rename Azimuth direction

bge <- bge %>% mutate(
  TA = paste(Tower,Antenna,sep=':'),
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

#df.all <- df.all %>% 
#  filter(grepl(pattern = '(^L[0-9]+-1$)|(^t[0-9]+-[0-9]+$)', x = L_ID)) %>% 
    # exclude L0 or L*-2 type
#  mutate_at(c('L_E','L_N','T_E','T_N'), floor) 
  
bge.median  <- bge %>% group_by(L_ID, Tower, Antenna, TA, 
    L_N, L_E, T_N, T_E, EstAntAz, CalcAzimuth, Distance, logDistance,
    AngleD, Angle, EstAzimuth, frontlobe, backlobe, anglef, angleb) %>%
    summarize(Power = median(Power)) %>% ungroup()

# omit all test points 

bge <- as.data.frame(bge)
bge.median <- as.data.frame(bge.median)

test.point <- c('L10-1', 'L13-1', 'L7-1', 't1-5', 't1-8',
  't2-4', 't2-8', 't3-5', 't3-9', 't4-9')

d.train <- subset(bge.median, !(L_ID %in% test.point)) 
# omit all test points from training data

# fit lmer model to training data

bge.mm <- lmer(Power ~ logDistance + anglef + angleb 
  + (1|Tower:Antenna) + (0+anglef|Tower:Antenna) 
  + (0+angleb|Tower:Antenna), subset = (Distance < 170),
  data=d.train)

# now estimate locations of test points


# bucket to save # obs, predicted location, and VC matrix
estloc <- matrix(NA, nrow=length(unique(test.point)), 
  ncol=8)

ntest <- length(unique(test.point))
for (i in 1:ntest) { 
  cat(i, ' ');
  
  bit <- subset(bge.median, L_ID == test.point[i])

# don't try to fit any location with only one tower
  if(length(unique(bit$TA)) == 1) {next}
  
  prediction <- optim(par = c(mean(bit$T_E), mean(bit$T_N)), 
    fn = ss.xy, hessian = T, 
    method = 'BFGS',
    ds=bit, model=bge.mm, trunc=54, s=5)
  
# VC matrix when using SS as the criterion.  temp$value is error SS  
  variance_matrix <- solve(prediction$hessian)*2*prediction$value/(n-2) 
    
     
# save results      
    estloc[i,] <- c(n=dim(test.distance)[1], prediction$par[1:2], 
      variance_matrix[c(1,2,4)], 
      bit$L_E[1], bit$L_N[1])
    }

dimnames(estloc)[[2]] <- c('Nobs','EstX','EstY', 
    'VarX','Cov','VarY', 'TrueX','TrueY')
 
   
# summarize results

delta <- sqrt((estloc[,2] - estloc[,7])^2 +
    (estloc[,3] - estloc[,8])^2)

median(delta)
range(delta)

nests <- dim(estloc)[1]
CIarea <- rep(NA, nests)

for (i in 1:nests) {
  thisest <- estloc[i,]
  m <- matrix(0, nrow=2, ncol=2)
  m[1,1] <- thisest[4]
  m[1,2] <- m[2,1] <- thisest[5]
  m[2,2] <- thisest[6]
  n <- thisest[1]
  CIarea[i] <- pi*sqrt(det(m))*(2*(n-1)*qf(0.95, 2, n-2)/(n-2))
}
effRadius <- sqrt(CIarea/pi)

# effective radius
range(effRadius)
median(effRadius)


