start.time <- Sys.time()
## generating a portfolio of VAs
library(vamc)
set.seed(123)
# n is the number of records
n <- 500
VAPort <- genPortInception(birthDayRng = c("1950-01-01", "1980-01-01"),
                           issueRng = c("2001-08-01", "2014-01-01"), matRng = c(15, 30),
                           acctValueRng = c(10000, 500000), femPct = 0, fundFee = c(30, 50, 60, 80, 10, 38, 45, 55, 57, 46), baseFee = 200, prodPct = rep(1/2, 2),
                           prodType = c("DBRP"), riderFee = c(90), rollUpRate = rep(5), withdrawalRate = 0, numPolicy = n) 
# extracting the numeric columns; maturity, age, premium
Maty <- list()
for (i in 1:n) {
  
  yy <- as.numeric(difftime(as.Date(VAPort[i,6]), as.Date(VAPort[i,5]), unit="weeks"))/52.2525
  Maty[[i]] <- round(yy)
  
}
Matyears <- matrix(unlist(Maty) , ncol = 1)

age <- list()
for (i in 1:n) {
  
  yy <- as.numeric(difftime(as.Date(VAPort[i,5]), as.Date(VAPort[i,7]), unit="weeks"))/52.2525
  age[[i]] <- round(yy)
}
Age <- matrix(unlist(age) , ncol = 1)

library(dplyr)
Fund <- select(VAPort, fundValue1, fundValue2, fundValue3, fundValue4, fundValue5, fundValue6, fundValue7, fundValue8, fundValue9, fundValue10 )
Premium <- matrix(rowSums(Fund), ncol = 1)
# combining three attributes age, premium, maturity in one matrix for clustering
clusterdata <- cbind(Age, Premium, Matyears)
x <- scale(clusterdata)

## finding the optimal number of clusters and clustering
library(stats)
library(factoextra)
fviz_nbclust(x, kmeans, method = "wss") + geom_vline(xintercept = 6 , linetype = 2)
# k is the number of cclusters
k <- 6
mymod <- kmeans(x=x, centers=k)

# Illusteration of clusters
#fviz_cluster(mymod, data = x,
#            ellipse.type = "euclid", # Concentration ellipse
#            star.plot = TRUE, # Add segments from centroids to items
#            repel = FALSE, # Avoid label overplotting (slow)
#            ggtheme = theme_minimal()
#)
x <- cbind(x,1:nrow(x)) 

#assigning each cluster to secific matrix c[[i]]
# the range of i is the number of clusters
cl = list()
for (i in 1:k) {
  cl[[i]] <- data.frame(x[mymod$cluster==i,])
}

library(dplyr)
y <- list()
disttt <- list()
# finding the distance of each record in each cluster from the center of the cluster
# range of r is the number of clusters
for (r in 1:k)
{
  # calculating the Euclidean distance of each column of each cluster from its cluster center
  # the range of j is the number of attributes or columns
  for(j in 1:3){
    y[[j]] <- (cl[[r]][,j] - mymod$centers[r,j])^2
  }
  
  #disttt[[r]] is the list of distances for the each record of each cluster r  
  disttt[[r]] <- sqrt(Reduce("+" , y))
  # attaching the calculated distances to the cluster matrix  
  cl[[r]] <- cl[[r]] %>% mutate(disttt[[r]])
}
#finding the minimum distance to find the representatives in the data set
result <- list()
for (j in 1:k) 
{
  result[[j]] <- cl[[j]][which.min(cl[[j]]$`disttt[[r]]`), ]
}

## Monte Carlo simulation and pricing the representatives
indexScen <- genIndexScen(mCov, 100, 360, indexNames, 1 / 12, cForwardCurve, 1)
indexScen[1, 1:5, ]
fundScen <- genFundScen(fundMap, indexScen)
fundScen[1, 1:5, ]

exPolicy <- list()
clvalue <- list()
for (i in 1:k) {
  pr <- result[[i]][,4]
  exPolicy[[i]] <- VAPort[pr, ]
  clvalue[[i]] <- valuateOnePolicy(exPolicy[[i]], mortTable, fundScen[1, , ], 1 / 12, cForwardCurve)
}

## "kriging"
Sum <- list()
Sum2 <- list()
Sum3 <- list()
# j and r are the number of clusters
#i is the number of attributes
for (j in 1:k) {
  for (r in 1:k) {
    for (i in 1:3) {
      # Euclidean distance for all clusters' repsetentatives      
      Sum[[i]]  <- (result[[j]][,i] - result[[r]][,i])^2
    }
    Sum2[[r]] <- sqrt(Reduce("+", Sum))
  }
  Sum3[[j]] <- Sum2 
}
# converting the list to matrix
V <- matrix(unlist(Sum3), ncol = k)

lum <- list()
lum2 <- list()
lum3 <- list()
# t is the number of records in the data set
# r is the number of clusters
# i is the number of attributes
for (t in 1:n) {
  for (r in 1:k)  {
    for (i in 1:3) {
      # Euclidean distance for each record from all r clusters' representatives     
      lum[[i]] <- (x[t,i] - result[[r]][,i])^2
    }
    lum2[[r]] <- sqrt(Reduce("+", lum))
  }
  lum3[[t]] <- lum2 
}
# converting list to matrix
D <- matrix(unlist(lum3), ncol = n)

# alpha and bet are the parameters difined in the article
beta <- quantile(V, probs = 0.95)
alpha <- 0
# final matrices of distances
Vjr <- alpha + exp((-3/beta)*V)
vc <- matrix(rep(1, k) , ncol = 1)
vr <- c(rep(1,k),0)
Vfin <- rbind(cbind(Vjr, vc),(vr))

Dtr <- alpha + exp((-3/beta)*D)
Dtotal <- rbind(Dtr, matrix(rep(1,n) , nrow = 1))
# calculating the matrix of weights for each record
W <- list()
for (i in 1:n) {
  W[[i]] <- solve(Vfin , Dtotal[,i])
  
}
# Wfin is the final matrix of weights for each record
Wfin <- matrix(unlist(W) , ncol = n)
Wxi <- Wfin[-(k+1),]
# Yr is the matrix of prices claculated by MC for representatives
Y <- matrix(unlist(clvalue) , nrow = 2)
Yr <- matrix(Y[2,] , ncol = 1)
# polvalue is the list of values for each record
polvalue <- list()
for (i in 1:n) {
  polvalue[[i]] <- t(Wxi[,i]) %*% Yr
  
}
## final portfolio value which is the sum of all values of records
PortfolioValue <- Reduce("+" , polvalue)

## Second Part: calculating weights for each cluster
# Dp is the distance of each cluster
Dp <- rbind(matrix(apply(Dtr, 1, sum) , ncol = 1) , n)
# Wc is the matrix of total weights of each cluster
Wc <- solve(Vfin, Dp)
# Portfolio value is the product of Wc and the price matrix of representatives
Wp <- matrix(Wc[-(k+1),], ncol = 1)

##final Portfolio value
PORTFOLIOVALUE <- t(Wp) %*% Yr
PORTFOLIOVALUE
end.time <- Sys.time()
time.taken <- round(end.time - start.time , 2)
time.taken