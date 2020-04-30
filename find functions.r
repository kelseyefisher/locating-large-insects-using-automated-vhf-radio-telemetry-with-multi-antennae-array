# functions used in find locations.r

# rewritten to use June 19 data structure

truncPower <- function(power, s, trunc=54) {
  # adjust power for truncation
  z <- (trunc-power)/s
  power + s*dnorm(z)/pnorm(z, lower=F)
  }

# tower / antenna information 
#   and an angle model with separate forward and backward lobes

plot.tower <- function(ds, ydata=NULL, label=F, 
  round=1, expansion = 1.1, addLoc=T, labelLoc=NULL, center=T) {
  # plot tower locations, antenna directions and (optionall) a vector response
  # antenna directions are clockwise degrees from N
  #   compute angles are counter clockwise degrees from E
  # center = F => use northing and easting;
  #   =T center plot at (0,0)
  
  x <- ds$T_E - mean(ds$T_E)
  y <- ds$T_N - mean(ds$T_N)
  xlim <- range(x)*expansion
  ylim <- range(y)*expansion
  
  if (!center) {
    x <- x + mean(ds$T_E)
    xlim <- xlim + mean(ds$T_E)
    y <- y + mean(ds$T_N)
    ylim <- ylim + mean(ds$T_N)
    }
  
  plot(x,y, xlim=xlim, ylim=ylim, type='n', 
    xlab='Easting', ylab='Northing', asp=1)
  
  if (label) {text(x, y, substring(ds$Tower, 2,2)) }
    else {points(x, y, pch='+', col=3, cex=1.2) }
  
  angle <- pi/2 - 2*pi*ds$EstAzimuth/360
  
  bitx <- 10*cos(angle)
  bity <- 10*sin(angle)
  
  segments(x, y, x + bitx, y + bity)
  if (!is.null(ydata)) {
    if (class(ydata) =='character') {
       text(x + 1.5*bitx, y + 1.5*bity, ydata)
       }
    else { 
      text(x + 1.5*bitx, y + 1.5*bity, round(ydata, round) )
      }
  }
  
  
  if (addLoc) {
    points(ds$L_E - mean(ds$T_E), ds$L_N - mean(ds$T_N), pch=19, col=4)
  }
  if (!is.null(labelLoc)) {
    text(ds$L_E - mean(ds$T_E), ds$L_N - mean(ds$T_N), 
      labelLoc)
  }

  invisible(NULL)
  }

merge.raneff <- function(ds, fit, inclTower=T, inclTA = T) {
  # merge random effects for Tower and/or Tower:Antenna from fit into the data set ds
  # uses dplyer functions
  
  ds$TAkey <- paste(ds$Tower, ds$Antenna, sep=':')
  # matches Tower:Antenna code in ranef(fit)
  
  # setup random effects vectors
  if (inclTower) {
    temp <- ranef(fit)$Tower
    rT <- data.frame(Tower=row.names(temp), ranTower=temp[,1], stringsAsFactors = F)
    ds <- left_join(ds, rT, by='Tower')
    }
  
  if (inclTA) {
    temp <- ranef(fit)$'Tower:Antenna'
    rTA <- data.frame(TAkey=row.names(temp), ranTA=temp[,1], stringsAsFactors = F)
    ds <- left_join(ds, rTA, by='TAkey')
    }
  ds
}

power <- function(xy, ds, model, trunc=54, s=5) {
  # predict power at each antenna for a transmitter at xy
  #   using expected truncated power with specified sd
  # parameter vector is (x, y)
  # ds: data set with tower locations, tower and antenna
  #   information for each observation 
  # model: lmer model to predict power
  # trunc: NULL when no adjustment for truncation
  
  x <- as.numeric(xy[1])
  y <- as.numeric(xy[2])

  theta <- atan2(x-ds$T_E, y-ds$T_N) - 2*pi*ds$EstAzimuth/360
    # angle in radians from proposed location to antenna boresight
 
  theta <- ((theta + pi) %% (2*pi)) - pi
    # centered at 0 down boresight, so -pi/2 is W, pi/2 is E
  AngleDs <- theta*360/(2*pi)

  frontlobe = (AngleDs > -90) & (AngleDs < 90)
  backlobe = 1-frontlobe
  anglef = frontlobe*((1-cos(theta))^(0.95)-1)
  angleb = backlobe*((1+cos(theta))^(0.95)-1)
  
  d <- sqrt(7^2 + (ds$T_E-x)^2 + (ds$T_N-y)^2)
    # distance, including 7m vertical
  
  predloc <- data.frame(Tower=ds$Tower, Antenna=ds$Antenna,
    logDistance=log(d), 
    anglef=anglef, angleb=angleb)
  
  estPower <- predict(model, newdata=predloc)
    # power from linear model 
  
  if (!is.null(trunc)) {
    # adjust for truncation, if requested
    
    alpha <- (trunc	- estPower)/s
    # numerical issues when alpha large (0/0), set tadj to 10*s
    tadj <- ifelse(alpha > 6, 
      10*s,
      s*(dnorm(alpha)/(1-pnorm(alpha))) )
    # numerical issues and when alpha small (< -10), set tadj to 0
    tadj <- ifelse(alpha < -6, 
      0,
      tadj )
      }
	 else{ tadj <- 0}
     # and set adjustment to 0 when not including truncation
  
  estPower + tadj
  }

lnl.xy4 <- function(xys, betaAngle, betaDist, ds, 
  inclTower=F, inclTA=F,
  inclDist=T, inclAngle=T, inclTrunc=T, trunc=54,
  debug=F) {
  
  sigma <- exp(xys[3])
  nobs <- dim(ds)[1]
  if (sigma > 35) {
    nlnl <- 100*nobs
    }
   # if sigma too big, return large neg lnl
  else {
  pred <- power(xys, betaAngle, betaDist, ds, 
    inclTower=inclTower, inclTA=inclTA, 
    inclDist=inclDist, inclAngle=inclAngle, 
    inclTrunc=inclTrunc, trunc=trunc
    )
  
  ss <- sum((pred - ds$Power)^2)
  
  nlnl <- 0.5*ss/sigma^2 + nobs*log(sigma)
  if (debug) {
    print(c(xys, ss, nlnl))
    }
  }
  nlnl
  }
  
ss.xy <- function(xy, ds, model, trunc=54, s=5) {
  # return SS fitting power in ds using model
  #   to predict at location in xy
  
  pred <- power(xy, ds, model, trunc, s)
  sum((pred - ds$Power)^2)
  
  }
