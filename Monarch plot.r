# Monarch plot.r

# plot the estimated locations and calculate delta and ciArea

# pdf(file='M100b.pdf', height=6, width=6)

par(mfrow=c(1,1), mar=c(3,3,0,0)+0.5, mgp=c(2,0.8,0))

plot.tower(Mmove1Median, label=F, center=F, 
  expansion=1.02)
points(Mtrue$L_E, Mtrue$L_N, pch=19, col=2)
points(estlocMedian[,1], estlocMedian[,2], 
  pch=19, col=c(1,4,4,4))
displaceX <- c(5, 5, 7, -20, -22)
displaceY <- c(4, 2, 0, -3, 0)
adj=c(0, 0,0,1, 0)
text(displaceX+c(estlocMedian[,1], Mtrue$L_E), 
  displaceY+c(estlocMedian[,2], Mtrue$L_N), 
  c('Observed', 'power+9.6', 'power+5', '+T3', 'True'),
  cex=0.7, adj=adj
  )


# add CI ellipse for obs location (row 1 of estlocMedian)
# ellipse is larger than entire area between towers
#   only visible if expansion in plot.tower ca 3

m <- matrix(NA, nrow=2, ncol=2)
m[1,1] <- estlocMedian[1,3]
m[1,2] <- m[2,1] <- estlocMedian[1,4]
m[2,2] <- estlocMedian[1,5]
n <- estlocMedian[1,6]
lines(ellipse(
  m, centre=estlocMedian[1,1:2],
  t = sqrt(2*(n-1)*qf(0.95, 2, n-2)/(n-2)) )
  )

# dev.off()

delta <- sqrt( (Mtrue$L_E - estlocMedian[,1])^2
  + (Mtrue$L_N - estlocMedian[,2])^2)
CIarea <- effRadius <- rep(NA, dim(estlocMedian)[1])
for (i in 1:dim(estlocMedian)[1]) {
  m <- matrix(NA, nrow=2, ncol=2)
  m[1,1] <- estlocMedian[i,3]
  m[1,2] <- m[2,1] <- estlocMedian[i,4]
  m[2,2] <- estlocMedian[i,5]
  n <- estlocMedian[i,6]
  CIarea[i] <- pi*sqrt(det(m))*(2*(n-1)*qf(0.95, 2, n-2)/(n-2))
  effRadius[i] <- sqrt(CIarea[i]/pi)
}
delta
effRadius
