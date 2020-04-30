# functions.r
#   Gang's functions - used in Prediction and plot.r

cal.az <- function(t_e, t_n, l_e, l_n,az.a = NULL){
  location <- c(l_e - t_e, l_n - t_n) 
  hypotenuse <- sqrt(location[1]^2 + location[2]^2)
  if (location[2] > 0) {angle <- acos(location[1]/hypotenuse) * 360/(2*pi)}
  else {angle <- 360 - (acos(location[1]/hypotenuse) * 360/(2*pi)) }
  angle <- 360 - angle + 90
  ifelse(angle > 360, angle <- angle - 360, angle <- angle)
  
  if (!is.null(az.a)){
    angle <- angle - az.a
    ifelse(angle < 0, angle <- angle + 360, angle <- angle)
    angle <- angle*2*pi/360
  }
  return(angle)
}


log.likelihood <- function (location){
  
  x <- location[1]; y <- location[2]
  
  test.distance <<- test.distance %>% mutate(
    Distance = sqrt((x - T_E)^2 + (y - T_N)^2),
    log_Dist = log(Distance),
    CalcAzimuth = Vectorize(cal.az)(T_E,T_N,x,y), 
    Az = ifelse(CalcAzimuth - EstAntAz < 0, CalcAzimuth - EstAntAz + 360, CalcAzimuth - EstAntAz),
    Az.2 = ifelse(Az > 180, Az - 360, Az),
    Angle_plus_mask = ifelse(cos(Az*2*pi/360) > 0, (1-cos(Az*2*pi/360))^Angle_mask_power-1, 0),
    Angle_minus_mask = ifelse(cos(Az*2*pi/360) < 0, (1+cos(Az*2*pi/360))^Angle_mask_power-1, 0))
  
  test.distance$predict <<- predict(model2.2,test.distance, allow.new.levels=T )
  
  test.distance <<- test.distance %>% mutate(OE = truncnorm::etruncnorm(a = 50, mean = predict, sd = sigma(model2.2)), 
    #ll = log(dnorm(Power, predict, sigma(model2.2))* (1-pnorm(50,predict, sigma(model2.2)))^(-1)),
    #ll= log(truncnorm::dtruncnorm(Power, a=50, b=Inf, mean = predict, sd = sigma(model2.2))), # same as above
    ll = log(dnorm(Power, predict, sigma(model2.2)))) 
  
  
  ss <- test.distance %>% transmute((Power - OE)^2) %>% sum()
  ll <- test.distance %>% pull(ll) %>% sum() %>% subtract()
  
  return(ll)  #  return(ss)  # here could change of ML method or LSS method
}

g2 <- function(theta){
  0.4*cos(abs(theta))+0.5*sqrt(4*0.4^2*(cos(abs(theta)))^2-4*(0.4^2-0.6^2))
}
#cut <- 90
g2_linear <- function(theta,theta0=cut){
  theta0 <- abs(theta0*pi/180)
  a <- (0.2-g2(theta0))/(pi-theta0)
  b <- 0.2-a*pi
  a * abs(theta) + b
}
#plot(g2, -pi, pi, col = "red")
#plot(g2_linear, -pi, pi, add=TRUE, col = "blue")

