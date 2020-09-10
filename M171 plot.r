# plots and summaries for M171 analysis and supporting material

# summary statistics

cat('Summary of distances btwn obs and estimated, obs power\n')
print( summary(alldist0) )
cat('Summary of distances btwn obs and estimated, power+22\n')
print( summary(alldist22) )

# distances for locations 19 and 59
cat('Distance btwn obs and estimated: obs power\n')
print( alldist0[c(19, 59)] )
cat('Distance btwn obs and estimated: power+22\n')
print( alldist22[c(19, 59)] )

# same computations for power as reported, uses allLocsB0
#   which is not computed here

# calculate effective radius of the CI for location 19
#  for power as reported, 
#  not run because allLocsB0 not computed here

allLoci <- allLocsB0[19,]
vc <- matrix(allLoci[c('Var E','Cov','Cov','Var N')], nrow=2)
n <- allLoci['Nobs']
  # F quantile for 95% confidence ellipse
q <- (n - 1) * 2 * qf(p = 0.95, df1 = 2, df2 = (n-2))/(n-2) 

cat('CI eff radius for loc 19, obs power\n')
print( sqrt(sqrt(det(vc*q)) ) )

# and for power with 22 added to compensate for veg interference
allLoci <- allLocsB22[19,]
vc <- matrix(allLoci[c('Var E','Cov','Cov','Var N')], nrow=2)
n <- allLoci['Nobs']
q <- (n - 1) * 2 * qf(p = 0.95, df1 = 2, df2 = (n-2))/(n-2) 
  # F quantile for 95% confidence ellipse

# effective radius
cat('CI eff radius for loc 19, power + 22\n')
print(sqrt(sqrt(det(vc*q)) ) )

# and repeat for location 59, first using obs power
allLoci <- allLocsB0[59,]
vc <- matrix(allLoci[c('Var E','Cov','Cov','Var N')], nrow=2)
n <- allLoci['Nobs']
q <- (n - 1) * 2 * qf(p = 0.95, df1 = 2, df2 = (n-2))/(n-2) 
  # F quantile for 95% confidence ellipse

cat('CI eff radius for loc 59, obs power\n')
print( sqrt(sqrt(det(vc*q)) ) )

# then using power with 22 added
allLoci <- allLocsB22[59,]
vc <- matrix(allLoci[c('Var E','Cov','Cov','Var N')], nrow=2)
n <- allLoci['Nobs']

# F quantile for 95% confidence ellipse
q <- (n - 1) * 2 * qf(p = 0.95, df1 = 2, df2 = (n-2))/(n-2) 

# effective radius
cat('CI eff radius for loc 59, power+22\n')
print(sqrt(sqrt(det(vc*q)) ) )

# plot of estimated locations

# this is the final version going into the paper
# colors identifying the two clusters of locations
clust171.col <- c(4,2)[1+(allLocsB22[,'est_E'] < 525600)]
true171.col <- c(4,2)[1+(trueloc171$Easting < 525645)]

# pdf(file='M171tower4.pdf', height=4, width=7)
par(mfrow=c(1,2), mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))
plot.tower(M171tower, bar=0, center=F, label=F,
  expansion=1.7)
points(allLocsB0[,'est_E'], allLocsB0[,'est_N'],
  pch=19, col=clust171.col)
points(trueloc171$Easting, trueloc171$Northing,
  col=true171.col)
add.ellipse(allLocsB0[19,], ellipse.col=4)
add.ellipse(allLocsB0[59,], ellipse.col=2)

legend('bottom', 'Using power as recorded', bty='n')

plot.tower(M171tower, bar=0, center=F, label=F,
  expansion = 1.7)
points(allLocsB22[,'est_E'], allLocsB22[,'est_N'],
  pch=19, col=clust171.col)
points(trueloc171$Easting, trueloc171$Northing, 
  col=true171.col)
add.ellipse(allLocsB22[19,], ellipse.col=4)
add.ellipse(allLocsB22[59,], ellipse.col=2)

legend('bottom', 'Using power+22', bty='n')

#dev.off()


