# data processing to set up variables for front- and back-lobe power

joint <- 170
Angle_mask_power <- 0.95
cut <- 90

df.all <- df.all %>% 
  filter(grepl(pattern = '(^L[0-9]+-1$)|(^t[0-9]+-[0-9]+$)', x = L_ID)) %>% 
    # exclude L0 or L*-2 type
  mutate_at(c('L_E','L_N','T_E','T_N'), floor) 
  
df.all <- df.all %>% 
  select(c(L_ID, Tower, Antenna, Time, Power, L_N, L_E, T_N, T_E, EstAntAz)) %>%
  mutate(Tower = gsub("[tT]", "", Tower),
           Antenna = as.character(Antenna),
           L_ID = as.character(L_ID),
           EstAntAz = EstAntAz + 0.7,
           CalcAzimuth = Vectorize(cal.az)(T_E,T_N,L_E,L_N)) %>% 
  mutate(Distance = sqrt((T_N-L_N)^2 + (T_E-L_E)^2),
           log_Dist = log(Distance),
           DS_1 = ifelse(Distance <= joint, Distance - joint, 0),
           DS_2 = ifelse(Distance <= joint, 0, Distance - joint),
           Az = ifelse(CalcAzimuth - EstAntAz < 0, CalcAzimuth - EstAntAz + 360, CalcAzimuth - EstAntAz),
           Az.2 = ifelse(Az > 180, Az - 360, Az),
           Az.2.r = Az.2*2*pi/360,
           Angle_plus_mask = ifelse(cos(Az*2*pi/360) > 0, (1-cos(Az*2*pi/360))^Angle_mask_power-1, 0),
           Angle_minus_mask = ifelse(cos(Az*2*pi/360) < 0, (1+cos(Az*2*pi/360))^Angle_mask_power-1, 0),
           AorR = ifelse(grepl("[tT]", L_ID), 'A', 'R'))
  
df.all <- df.all %>% 
  full_join(data.frame(Tower = rep(1:4 %>% as.character(), each = 4), 
     Antenna = rep(1:4 %>% as.character(),4), T_A = as.character(1:16)), 
    by = c('Tower','Antenna')) %>% 
  mutate(T_A = as.character(T_A)) %>% 
  dplyr::group_by(L_ID, T_A) %>% 
  mutate(sd = sd(Power)) %>% ungroup() %>% as.data.frame() 
  
df.all$EstAzimuth <- df.all$EstAntAz
  #  make copy of Antenna azimuth using my variable name
  
pool_sd <- df.all %>% dplyr::group_by(L_ID, T_A) %>% 
    mutate(n = n()) %>% ungroup() %>% 
    filter(n > 1) %>% mutate(n=NULL) %>% 
    group_by(L_ID) %>% 
    mutate(pool_sd = ANOVAreplication::pooled.sd(data.frame(Power, T_A))) %>% 
    select(L_ID,pool_sd) %>% unique() %>% as.data.frame()
  
df.all <- df.all %>% full_join(pool_sd, by = 'L_ID') %>% as.data.frame()


df.all <- df.all %>% mutate(g2 = g2(Az.2.r))
  
df.all <- df.all %>% mutate(g2_inside = ifelse(between(Az.2, -abs(cut),abs(cut)), g2-g2(abs(cut)*pi/180), 0),
                     g2_outside = ifelse(between(Az.2, -abs(cut),abs(cut)), 0, g2-g2(abs(cut)*pi/180)),
                     g2_outside_linear = ifelse(between(Az.2, -abs(cut),abs(cut)), 0, g2_linear(Az.2.r, theta0 = cut) - g2(abs(cut)*pi/180)))
  
  df.all <- df.all %>% mutate(g3 = g2_linear(theta = Az.2.r, theta0 = cut))

  median.all <- df.all %>% group_by(L_ID, Tower, Antenna, 
    L_N, L_E, T_N, T_E, T_A, EstAntAz, CalcAzimuth, Distance, log_Dist,
    Angle_plus_mask, Angle_minus_mask, Az.2, Az.2.r, pool_sd, 
    g2, g2_inside, g2_outside, g2_outside_linear, g3) %>%
    summarize(Power = median(Power))
  
df.all <- median.all %>% ungroup()
  # use median power at a location
  
median.all$EstAzimuth <- median.all$EstAntAz
  # and make a copy of angle variable for plot.tower()

data.070718.all <- as.data.frame(median.all)
  # use median power values 

data.070718.all$Date <- '070718'; data.070718.all$Field <- 'BGE'

# omit all test points 

test.point <- c('L10-1', 'L13-1', 'L7-1', 't1-5', 't1-8',
  't2-4', 't2-8', 't3-5', 't3-9', 't4-9')

d.train <- subset(data.070718.all,
    Distance > 0 & Distance < 190 & !(L_ID %in% test.point)) 
# omit all test points from training data

# fit lmer model to training data

initial.x <- mean(unique(d.train$T_E), na.rm = T)
initial.y <- mean(unique(d.train$T_N), na.rm = T)

model2.2 <- lmer(Power ~ log_Dist + Angle_plus_mask + Angle_minus_mask + 
                     (0+Angle_plus_mask|T_A) + (0+Angle_minus_mask|T_A) + (1|T_A),
					 subset=Distance < 170,
                   data = d.train, REML = T)

# now estimate locations of test points

# bucket to save # obs, predicted location, and VC matrix
estloc <- estnls <- matrix(NA, nrow=length(unique(test.point)), 
  ncol=8)

for (i in  seq_along(unique(test.point))){ 
  cat(i, ' ');
  
  test.distance <- data.070718.all %>% filter(L_ID == test.point[i])
  test.distance$OE <- vector(mode = "numeric", length = nrow(test.distance))
  test.distance$ll <- vector(mode = "numeric", length = nrow(test.distance))

# don't try to fit any location with only one tower
  if(length(unique(test.distance$T_A)) == 1) {next}
  
  prediction <- optim(par = c(initial.x, initial.y), 
    fn = log.likelihood, hessian = T, 
    method = 'L-BFGS-B')
  
  x <- prediction$par[1]
  y <- prediction$par[2]
  x.true <- test.distance$L_E[1]
  y.true <- test.distance$L_N[1]
  
  variance_matrix <- solve(prediction$hessian) 
    
     
# save results      
    estloc[i,1] <- dim(test.distance)[1]
    estloc[i,2] <- x
    estloc[i,3] <- y
    estloc[i,4] <- variance_matrix[1,1]
    estloc[i,5] <- variance_matrix[1,2]
    estloc[i,6] <- variance_matrix[2,2]
    estloc[i,7] <- x.true
    estloc[i,8] <- y.true
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


