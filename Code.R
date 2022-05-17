library(raster)
library(Matrix)
library(stocc)
library(stats)
library(boot)
library(spatstat)
library(fields)
library(RColorBrewer)
library(rasterVis)
library(truncnorm)
library(mvtnorm)
library(gridExtra)
library(fBasics)
library(coda)
library(beanplot)
library(mgcv)
library(gbm)
library(rgdal)

###Simulation Study
nsite=6400
nobs=2
nyear=16

#create cell layer for analyses
x <- raster(nrows=80, ncols=80, xmn=0, xmx=80, ymn=0,ymx=80, resolution=1, vals=1:nsite)
writeRaster(x, filename = "cell_sim.tif", "GTiff", overwrite=TRUE)
cell <- raster("cell_sim.tif")
plot(cell)


##simulate random positives for year 1
generateGaussianData <- function(n, center, sigma) {
  data = rmvnorm(n, mean = center, sigma = sigma)
  data = data.frame(data)
  names(data) = c("x", "y")
  # data = data %>% mutate(class=factor(label))
  data
}

n.pos = 5
center = c(40, 40)
sigma = matrix(c(2, 0, 0, 2), nrow = 2)
data1 = generateGaussianData(n.pos, center, sigma)

simgrid <- expand.grid(x=1:80,y=80:1)
rasValue=extract(x, data1)
pos.1 <- data.frame(simgrid,data= rep(0,times=nsite))
pos.1[rasValue,3] <- 1

##insert positives into matrix of truth
vec.matrix <- numeric(0)
true.presence <- matrix(vec.matrix, nrow = nsite, ncol = nyear)
true.presence[,1] <- pos.1[,3]

psi.pred <- matrix(vec.matrix, nrow = nsite, ncol = nyear)

##create spatial key
key <- data.frame(site=1:nsite, X=simgrid[,1]-0.5, Y=simgrid[,2]-0.5)
xy <- c(key$X, key$Y)
xy <- matrix(xy,  c(nsite, 2))


#function to ID neighbors in queen's configuration
neighborhood <- function(raster){
  nn <- matrix(,length(raster[]),8)
  for(i in 1:dim(nn)[1]){
    neigh <- adjacent(cell, i, directions=8, pairs=FALSE)
    fill <-length(neigh)
    fill.na <- rep(NA, 8-fill)
    nn[i,] <- c(neigh,fill.na)
  }
  nn
}

#run neighborhood function of study area grid
nn <- neighborhood(cell)

#effect sizes for state processes
alpha.est <- -3.7
neigh.occ <- .5
alpha.emer <- -11
alpha.pers <- 2.4
dist.emer <- -5.5
hab.cov.emer <- 0.5
hab.cov.pers <- -1.7

edge <- c(which(key[,2]> 78), which(key[,3]> 78), which(key[,2]<2), which(key[,3]<2))
edge <- sort(edge)
edge <- unique(edge)

# # ##Simulated RSR
rmvn <- function(n, mu = 0, V = matrix(1)) {
  p <- length(mu)
  if (any(is.na(match(dim(V), p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p)))
}
#
n <- nrow(simgrid)
distance <- as.matrix(dist(simgrid))
phi <- 0.3
habdat.pers.raw <- rmvn(1, rep(0, n), exp(-phi * distance))
habdat.emer.raw <- rmvn(1, rep(0, n), exp(-phi * distance))

habdat.pers <- ((habdat.pers.raw)-mean(habdat.pers.raw))/sd(habdat.pers.raw)
habdat.pers <- as.vector(habdat.pers)

habdat.emer <- ((habdat.emer.raw)-mean(habdat.emer.raw))/sd(habdat.emer.raw)
habdat.emer <- as.vector(habdat.emer)

habdat.emer[edge] <- min(habdat.emer)
habdat.pers[edge] <- min(habdat.pers)

##make colonization matrix

dist.std.pred <- matrix(vec.matrix, nrow = nsite, ncol = nyear) 
i.nbors.out <- matrix(vec.matrix, nrow = nsite, ncol = nyear)

#forecasting code
for(l in 1:(nyear-1)){
  
  #make colonization matrix
  true.presence.neigh <- data.frame(site=1:nsite, data=true.presence[ ,l])
  WhichNborsMat.pred <- true.presence.neigh$data[match(nn, true.presence.neigh$site)]
  WhichNborsMat.pred[is.na(WhichNborsMat.pred)] <- 0
  WhichNborsMat.pred <- matrix(WhichNborsMat.pred, nrow = nsite)
  
  # forecasted neighborhood colonization, n.neigh= number of infected neighbors, i.nbors= indicator if neighbors are infected (binary:0,1).
  n.neigh.pred <- rowSums(WhichNborsMat.pred);
  i.nbors.pred <- ifelse(n.neigh.pred > 0, 1, 0)
  i.nbors.out[,l] <- i.nbors.pred
  
  #sites of positives for each iteration
  positives.site <- true.presence.neigh[which(true.presence.neigh$data==1), ]
  allpoints.ref <- data.frame(site=1:nsite,x=key$X, y=key$Y)
  positives.x <- allpoints.ref$x[match(positives.site$site, allpoints.ref$site)]
  positives.y <- allpoints.ref$y[match(positives.site$site, allpoints.ref$site)]
  positives <- cbind(positives.x,positives.y)
  positives <- ppp(positives[,1], positives[,2], c(0, 80), c(0, 80))
  allpoints <- cbind(key$X, key$Y)
  allpoints <- ppp(allpoints[,1], allpoints[,2], c(0, 80), c(0, 80))
  distance <- nncross(allpoints,positives)
  distance <- distance$dist
  
  #standardize distance
  dist.std.pred[,l]=((distance)-mean(distance))/sd(distance)
  
  gamma <- inv.logit(alpha.est + neigh.occ * n.neigh.pred)  #established disease spread
  delta <- inv.logit(alpha.emer + dist.emer * dist.std.pred[ ,l] + hab.cov.emer * habdat.emer) #emergent disease spread
  phi <- inv.logit(alpha.pers + hab.cov.pers * habdat.pers)
  psi.pred[,l+1] <- (true.presence[ , l] * phi) + ((1-true.presence[ , l]) * i.nbors.pred * gamma) + ((1-true.presence[ , l]) * (1- i.nbors.pred) * delta)      #gamma established spread, delta emergent spread
  true.presence[ ,(l+1)] <- rbinom(length(psi.pred[,l+1]), size = 1, prob=psi.pred[,l+1])
  rm(positives.site, positives.x, positives.y, positives, true.presence.neigh, WhichNborsMat.pred, n.neigh.pred, i.nbors.pred, distance)
} #l

psi.pred[,1] <- psi.pred[,2]

###effect sizes for observation processes
alpha.det <- -1.5
weight.det <- 0.2
prev.det <- 0.5
mean.weight <- 0.4
sd.weight <- 1.8

#simulate prevalence trends over time
prev <- numeric(0)
for (l in 1:nyear){
  prev[l] <- 0.01*(1+.19)^l
  
}

prev.std <- ((prev)-mean(prev))/sd(prev)

prev.std.mat <- array(c(vec.matrix),dim = c(nsite,nyear))

for (l in 1:nyear){
  prev.std.mat[,l] <- rep(prev.std[l], times=nsite)
}

#simulate disease surveillance weights
weights <- rtruncnorm(1e6, a=0,b=Inf, mean=2,sd=5)


# #sample according to risk
multiplier.array <- array(c(vec.matrix),dim = c(nsite,nobs,nyear))

for(l in 1:nyear){
  
  multiplier.ar <- rep(psi.pred[,nyear], times=2)
  multiplier.array[,,l] <- matrix(multiplier.ar, nrow = nsite, ncol = 2)
  rm(multiplier.ar)
}


#increase probability that site is sampled beyond occurrence probability#
multiplier.array <- multiplier.array*10
multiplier.array <- ifelse(multiplier.array >= 1, 0.99, multiplier.array)

multiplier <- rbinom(length(multiplier.array), size = 1, prob=multiplier.array)
index <- which(multiplier==1)
index.total <- c(index)
index.total <- unique(index.total)
sample.weight <- sample(weights, length(index.total), replace=TRUE)
multiplier[index.total] <- sample.weight
w2 <- array(c(multiplier),dim = c(nsite,nobs,nyear))

#standardize
w2.std.y=(w2-mean.weight)/sd.weight

#create matrix to store data
vec.matrix <- numeric(0)
det.prob <- array(c(vec.matrix),dim = c(nsite,nobs,nyear))
eff.det.prob <- array(c(vec.matrix),dim = c(nsite,nobs,nyear))
y <- array(c(vec.matrix),dim = c(nsite,nobs,nyear))

for(l in 1:nyear){
  for(k in 1:nobs){
    det.prob[,k,l] <- inv.logit(alpha.det + prev.std.mat[,l] * prev.det + weight.det * w2.std.y[,k,l])
    eff.det.prob[,k,l] <- true.presence[,l] * det.prob[,k,l]
    y[,k,l] <- rbinom(n = nsite, size = 1, prob = eff.det.prob[,k,l])
    
  }
} 

##generate neighborhood data
pos <- matrix(vec.matrix, nrow = nsite, ncol = nyear) 
##Create a reference dataframe
for(l in 1:nyear){
  pos[,l] <- rowSums(y[,,1:l])
}

pos <- ifelse(pos[,] > 0, 1, 0)

WhichNborsMat.y <- array(vec.matrix, c(nsite, 8, nyear))
dist.std.y <- matrix(vec.matrix, nrow = nsite, ncol = nyear)
dist.y <- matrix(vec.matrix, nrow = nsite, ncol = nyear)

#forecasting code
for(l in 1:(nyear)){
  
  #make colonization matrix
  true.presence.neigh <- data.frame(site=1:nsite, data=pos[ ,l])
  WhichNborsMat.pred <- true.presence.neigh$data[match(nn, true.presence.neigh$site)]
  WhichNborsMat.pred[is.na(WhichNborsMat.pred)] <- 0
  WhichNborsMat.y[,,l] <- matrix(WhichNborsMat.pred, nrow = nsite)
  
  #sites of positives for each iteration
  positives.site <- true.presence.neigh[which(true.presence.neigh$data==1), ]
  allpoints.ref <- data.frame(site=1:nsite,x=key$X, y=key$Y)
  positives.x <- allpoints.ref$x[match(positives.site$site, allpoints.ref$site)]
  positives.y <- allpoints.ref$y[match(positives.site$site, allpoints.ref$site)]
  positives <- cbind(positives.x,positives.y)
  positives <- ppp(positives[,1], positives[,2], c(0, 80), c(0, 80))
  allpoints <- cbind(key$X, key$Y)
  allpoints <- ppp(allpoints[,1], allpoints[,2], c(0, 80), c(0, 80))
  distance <- nncross(allpoints,positives)
  distance <- distance$dist
  dist.y[,l] <- distance
  
  #standardize distance
  dist.std.y[,l]=((distance)-mean(distance))/sd(distance)
  rm(positives.site, positives.x, positives.y, positives, true.presence.neigh, distance)
  
}


#####Subset all for nyears######

nyear.sub <- 16
y.sub <- y[,,2:nyear.sub]
WhichNborsMat.y <- WhichNborsMat.y[,,2:nyear.sub]
dist.std.y <- dist.std.y[,2:nyear.sub]
w2.std.y <- w2.std.y[,,2:nyear.sub]
prev.std.mat <- prev.std.mat[,2:nyear.sub]


####Spatial Autocorrelation############################
#######################################################

occupancy.model <- model.matrix(~1,data=data.frame(rep(1,nsite)))
site=1:nsite

spatial.model=list(model="rsr", moran.cut=250)

Xz <-as.matrix(occupancy.model, site)
Xz.vec <- as.vector(Xz)
n.site=nrow(Xz)
Q <- icar.Q(xy, threshold= 1.42, rho = 1)
A <- diag(diag(Q)) - Q
P <- diag(n.site) - Xz %*% solve(crossprod(Xz), t(Xz))
Op <- (nrow(A)/sum(A)) * (P %*% (A %*% P))
e <- rARPACK::eigs(Op, as.integer(spatial.model$moran.cut))
K <- e$vectors
KtK <- diag(ncol(K))
Q.alpha <- as.matrix(t(K) %*% Q %*% K)
mu.a <- rep(0, nrow(Q.alpha))


#create list to load into jags
test1 <- list(y=y.sub, nsite=nrow(y.sub), nyear = dim(y.sub)[3], nrep = dim(y.sub)[2],
              WhichNborsMat = WhichNborsMat.y, dist=dist.std.y,
              K=K, Q.alpha= Q.alpha, mu.a=mu.a, Xz.vec=Xz.vec, hab.emer=habdat.emer,
              hab.pers=habdat.pers, W2=w2.std.y, prev=prev.std.mat)


rm(list = ls()[!ls() %in% c("test1")])

str(test1)
attach(test1)

# Specify model in BUGS language
sink("mod.txt")
cat("
      
      model {
      
      #Likelihoods
      
      # Observation model(yr1)
      for (i in 1:nsite){
      for (j in 1:nrep){
      muy[i,j,1] <- z[i,1] * p[i,j,1]
      logit(p[i,j,1]) <- alpha.det + weight.det * W2[i,j,1] + prev.det * prev[i,1]
      y[i,j,1] ~ dbern(muy[i,j,1])
      
      } #j
      
      #state model(yr1)
      
      z[i,1] ~ dbern(psi[i,1])
      logit(psi[i,1]) <- Xz.vec[i] %*% alpha.occ + RSR[i]
      } #i
      
      #observation for subsequent years
      
      for (k in 1:(nyear-1)){
      for (i in 1:nsite){
      for (j in 1:nrep){
      muy[i,j,(k+1)] <- z[i,(k+1)] * p[i,j,(k+1)]
      logit(p[i,j,(k+1)]) <- alpha.det + weight.det * W2[i,j,(k+1)] + prev.det * prev[i,(k+1)]   #add supplementary weight here
      y[i,j,(k+1)] ~ dbern(muy[i,j,(k+1)])
      } #j
      
      z.forjags[i,k] <- z[i,k]
      
      # #Neighborhood colonization, n.neigh= number of infected neighbors, i.nbors= indicator if neighbors are infected (binary:0,1).
      n.neigh[i, k] <- sum(WhichNborsMat[i, ,k]);
      i.nbors[i, k] <- ifelse(n.neigh[i, k] > 0, 1, 0)
      
      #Psi Mixture
      psi[i, (k+1)] <- z.forjags[i,k] * phi[i,k] + (1-z.forjags[i,k]) * i.nbors[i,k] * gamma[i,k] + (1-z.forjags[i,k]) * (1- i.nbors[i,k]) * delta[i,k]  #gamma established spread, delta emergent spread
      logit(gamma[i,k]) <- alpha.est + neigh.occ * n.neigh[i,k]  #established disease spread
      logit(delta[i,k]) <- alpha.emer + dist.emer * dist[i,k] + hab.cov.emer * hab.emer[i]  #emergent disease spread
      logit(phi[i,k]) <- alpha.pers + hab.cov.pers * hab.pers[i]
      z[i,(k+1)] ~ dbern(psi[i,(k+1)])
      } #i
      
      } #k
      
      for (k in 1:nyear){
      
      n.occ[k]<-sum(z[1:nsite,k])
      } #k
      
      
      #Priors
      #Priors
      alpha.occ ~ dunif(-10, 10)    #intercept for state model
      alpha.det ~ dunif(-10, 10)    #intercept for detection model
      weight.det ~ dunif(-10, 10)    #beta for effect of weighted probability on detection
      alpha.emer ~ dunif(-10, 10)   #intercept for long-distance spread
      dist.emer ~ dunif(-10, 10)    #beta for effect of distance from prior positives on state model
      alpha.est ~ dunif(-10, 10)    #intercept for localized spread
      neigh.occ ~ dunif(-10, 10)    #beta for effect of number of infected neighbors on state model
      alpha.pers ~ dunif(-10, 10)    #intercept for persistence
      hab.cov.emer ~ dunif(-10, 10)   #beta for effect of habitat covariate on long-distance spread
      hab.cov.pers ~ dunif(-10, 10)   #beta for effect of habitat covariate on persistence
      prev.det ~ dunif(-10,10)        #beta for effect of prevalence on detection

      #RSR Random spatial correlation
      sigma ~ dunif(0,100)
      tau <- pow(sigma, -2)
      alpha ~ dmnorm(mu.a, tau * Q.alpha)
      RSR <- K %*% alpha
      
      }
      
      ",fill = TRUE)
sink()


# Bundle data
win.data <- list(y = y, nsite = dim(y)[1], nrep = dim(y)[2], nyear = dim(y)[3], 
                 WhichNborsMat = WhichNborsMat, dist=dist, K=K, Q.alpha=Q.alpha, mu.a=mu.a, Xz.vec=Xz.vec, hab.emer=hab.emer, hab.pers=hab.pers, W2=W2, prev=prev)

rm(test1)

# Initial values
zst <- apply(y, c(1, 3), max)	# Observed occurrence as inits for z
inits <- function(){ list(z = zst)}

# Parameters monitored
params <- c( "alpha.occ", "dist.emer", "alpha.det","neigh.occ", "alpha.est", "alpha.emer", "alpha.pers", "n.occ", "psi", "hab.cov.emer", "hab.cov.pers", "prev.det", "weight.det")

# MCMC settings
ni <- 5000
nt <- 5
nb <- 50
nc <- 3

#Load the correct library
library("jagsUI")

# Call JAGS
out <- jags(win.data, inits, params, "mod.txt", n.chains = nc, n.thin = nt,
            n.iter = ni, n.burnin = nb)
