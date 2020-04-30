# plot BGE predictions and true locations
# assumes BGE prediction.r has been run 

# the important values are in the estloc data frame
 
# find bounding box for towers
 
xlim_low <- floor(min(c(d.train$T_E, d.train$L_E), na.rm= T)/100)*100
xlim_top <- ceiling(max(c(d.train$T_E, d.train$L_E), na.rm= T)/100)*100
ylim_low <- floor(min(c(d.train$T_N, d.train$L_N), na.rm= T)/100)*100
ylim_top <- ceiling(max(c(d.train$T_N, d.train$L_N), na.rm= T)/100)*100

# plot the towers and the antennae directions
par(mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))
plot.tower(median.all, addLoc=F, center=F)


# plot observed and estimated locations and connect with lines
estloc <- as.data.frame(estloc)
points(estloc$TrueX, estloc$TrueY,  col = "black")
points(estloc$EstX, estloc$EstY, col = "red")
segments(estloc$TrueX, estloc$TrueY, estloc$EstX, estloc$EstY, lty = 1, lwd = 1)

# add ellipses
npts <- dim(estloc)[1]
for (i in 1:npts) {
  thisloc <- estloc[i,]
  true <- thisloc[c('TrueX','TrueY')]
  xy <- thisloc[c('EstX','EstY')]
  diff <- t(as.matrix(true - xy))
  variance_matrix <- matrix(
    c(thisloc$VarX, thisloc$Cov, thisloc$Cov, thisloc$VarY), 
    ncol=2)
  Ivar <- solve(variance_matrix)
  if (variance_matrix [1,1] > 0 & variance_matrix[2,2] > 0){
      dir <- eigen(variance_matrix)[[2]][,1] # this define the direction of long axes of ellipse
      ag <- 360 - cal.az(0,0,dir[1], dir[2]) + 90
      
      n <- sum(!is.na(test.distance$Power))
      q <- (n - 1) * 2 * qf(p = 0.95, df1 = 2, df2 = (n-2))/(n-2) 
	  # old critical value of the scaled F distribution
      a.2 <- sqrt(eigen(variance_matrix)[[1]] * q)[1]
      b.2 <- sqrt(eigen(variance_matrix)[[1]] * q)[2]
      inorout <- (t(diff) %*% Ivar %*% diff) <= q
      ellipse.col <- ifelse(inorout,  "black", "red")
	  lines(ellipse(variance_matrix, centre=t(xy), t = sqrt(q) ),
  	    col=ellipse.col)
      }
   }
   
  

