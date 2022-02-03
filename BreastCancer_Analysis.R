##########################################
#  Analysis of Breast Cancer Dataset
#
#  Skills Used:                            
##########################################


############
# Read in the data
############
data = read.csv("Breast_Cancer.csv", header=T)
head(data)  



############
# center and scale the data
############
datanew = data[,-c(1,2)]  # get a dataframe with only the 30 variables that we are interested in
head(datanew)  

X = scale(datanew, scale=T) # center and scale the data
head(X)

D = dist(X)   # get the scaled distance matrix
head(D)
class(D) # for hierarchical clustering need class "dist"



###########
# Non-probabilistic k-means
###########
#### choose initial values  for centers of two clusters 
max(D)  # tells us the largest distance between rows
D1=as.matrix(D)  #Need this to find rows most distant.

rows.most.distant <- matrix(c(0, 0), ncol = 1) 
#### Note: the following nested for loop takes a really long time to run, so I commented it out
#### But you can uncomment it and run it to see that I get 213 & 153 as my most distant rows
# for(i in 1:569){
#   for(j in 1:569){
#     if(D1[i, j] == max(D1)  ) {
#       rows.most.distant <- matrix(c(i, j), ncol = 1)
#     }
#   }
# }
# rows.most.distant # rows 213 & 153 are the most distant. We set these as initial values. 


# Initial means for the 30 variables in cluster 1 
c1=X[213,] #initial center for cluster 1
c1
#Initial means for the 30 variables in cluster 2 
c2=X[153,]  #initial center for cluster 2 
c2
# Two index vectors (Z vectors): one  will tell us clusters assigned in last iteration
# Another index vector will tell cluster in the new iteration
# We make the indicator vectors as different as possible to start. 
pastIndicator=569:1   #initial value for z 
indicator=1:569     # past indicator will be compared with new indicator
### note: we initialize this way to get the algorithm started

###### We must iterate until pastIndicator=indicator
## While the two indicator vectors are different, keep going. 
while( sum(pastIndicator!=indicator)!= 0 ) 
{
  pastIndicator=indicator
  
  #distance to current cluster centers
  dc1 = colSums((t(X)-c1)^2)   #distance of each row from mean of cluster 1
  dc2 = colSums((t(X)-c2)^2)   #distance of each row from mean of cluster 2
  dMat = matrix(c(dc1,dc2),,2)
  
  #decide which cluster each point belongs to 
  indicator = max.col(-dMat)
  
  # update the cluster centers
  c1 = colMeans(X[indicator==1,])
  c2 = colMeans(X[indicator==2,])
}

##See the means of each cluster. 
c1   
c2   

## double check that Z_t-1 = Z(t) (convergence)
sum(indicator == pastIndicator) == nrow(datanew) # if we have reached convergence, this should return TRUE

# Use the last Z as indicator for the cluster 
cluster_allocation = indicator
head(cluster_allocation)
table(cluster_allocation)

###########
# Tables & Scatterplots for nonprobabilistic k-means clustering
###########

#### Extract the diagnosis of observations in our data
diagnosis = factor(data[,2])
head(diagnosis)

#### Get a table with which diagnosis falls in each cluster
#### We can also see if there is any misclassification
table(cluster_allocation, diagnosis)  
# Note that there is a pretty clear distinction that most of the benign diagnoses fall in cluster 2, and
# most of the malignant diagnoses fall in cluster 1.  
# However, there are 50 misclassifications in total
# So our classification using nonprobabilistic k-means was 91.2% accurate

#### Now we get a table of cluster percentage in each diagnosis
prop_c1_B <- 13/(344+13)
prop_c2_B <- 344/(344+13)
prop_c1_M <- 175/(175+37)
prop_c2_M <- 37/(175+37)
data.frame("cluster" = c(1,2), 'benign'=c(prop_c1_B,prop_c2_B),
           'malignant' = c(prop_c1_M,prop_c2_M))

#### Get a table of diagnosis percentage in each cluster
prop_B_c1 <- 13/(13+175)
prop_M_c1 <- 175/(13+175)
prop_B_c2 <- 344/(344+37)
prop_M_c2 <- 37/(344+37)
data.frame("cluster" = c(1,2), 'benign'=c(prop_B_c1,prop_B_c2),
           'malignant' = c(prop_M_c1,prop_M_c2))


## Get the scatterplots of each cluster plotted against 2 variables
#### mean radius & mean texture
plot(data[,3], data[,4],
     col=c("red","blue")[unclass(cluster_allocation)],pch=c(20,23)[unclass(diagnosis)],
     main="Nonprobabilistic k-means Clustering of Breast Cancer data, with mean radius vs mean texture",
     xlab="Mean radius (in mm)", ylab="Mean Texture")
legend("topleft",c("Benign","Malignant"),pch=c(20,23) )
legend("bottomright",c("cluster 1","cluster 2"),pch=c("R","B"),col=c("red","blue"))

#### mean perimeter & mean area  
plot(data[,5], data[,6],
     col=c("red","blue")[unclass(cluster_allocation)],pch=c(20,23)[unclass(diagnosis)],
     main="Nonprobabilistic k-means Clustering of Breast Cancer data, with mean perimeter vs mean area",
     xlab="Mean perimeter (in mm)", ylab="Mean area (in mm^2)")
legend("topleft",c("Benign","Malignant"),pch=c(20,23) )
legend("bottomright",c("cluster 1","cluster 2"),pch=c("R","B"),col=c("red","blue"))

#### mean smoothness & mean compactness  
plot(data[,7], data[,8],
     col=c("red","blue")[unclass(cluster_allocation)],pch=c(20,23)[unclass(diagnosis)],
     main="Nonprobabilistic k-means Clustering of Breast Cancer data, with mean smoothness vs mean compactness",
     xlab="Mean smoothness", ylab="Mean compactness")
legend("topleft",c("Benign","Malignant"),pch=c(20,23) )
legend("bottomright",c("cluster 1","cluster 2"),pch=c("R","B"),col=c("red","blue"))

#### mean concavity & mean concave points
plot(data[,9], data[,10],
     col=c("red","blue")[unclass(cluster_allocation)],pch=c(20,23)[unclass(diagnosis)],
     main="Nonprobabilistic k-means Clustering of Breast Cancer data, with mean concavity vs mean concave points",
     xlab="Mean concavity", ylab="Mean concave points")
legend("topleft",c("Benign","Malignant"),pch=c(20,23) )
legend("bottomright",c("cluster 1","cluster 2"),pch=c("R","B"),col=c("red","blue"))

#### mean symmetry & mean fractal dimension
plot(data[,11], data[,12],
     col=c("red","blue")[unclass(cluster_allocation)],pch=c(20,23)[unclass(diagnosis)],
     main="Nonprobabilistic k-means Clustering of Breast Cancer data, with mean symmetry vs mean fractal dimension",
     xlab="Mean symmetry", ylab="Mean fractal dimension")
legend("topleft",c("Benign","Malignant"),pch=c(20,23) )
legend("bottomright",c("cluster 1","cluster 2"),pch=c("R","B"),col=c("red","blue"))

##########
# Summary statistics for each cluster
#########
#### summary for cluster 1:  
mean_1 = c1 # mean of cluster 1, given by nonprobabilistic k-means algorithm
mean_1 

median_1 <- apply(datanew[cluster_allocation == 1, ], 2, median) # median of cluster 1
median_1  

var_covar_1 <- (1 / nrow(datanew[cluster_allocation == 1,])) *
  t(datanew[cluster_allocation == 1, ]) %*% scale(datanew[cluster_allocation == 1, ], scale=F)  # variance-covariance matrix of cluster 1
# we know the diagonal elements of var-covar matrix provide the variance of individual objects
var_1 <- diag(var_covar_1) # variance of cluster 1  
var_1

sd_1 <- sqrt(var_1) # standard deviation of cluster 1
sd_1

#### summary for cluster 2:
mean_2 = c2 # mean of cluster 1, given by nonprobabilistic k-means algorithm
mean_2 

median_2 <- apply(datanew[cluster_allocation == 2, ], 2, median) # median of cluster 1
median_2  

var_covar_2 <- (1 / nrow(datanew[cluster_allocation == 2,])) *
  t(datanew[cluster_allocation == 2, ]) %*% scale(datanew[cluster_allocation == 2, ], scale=F)  # variance-covariance matrix of cluster 1
# we know the diagonal elements of var-covar matrix provide the variance of individual objects
var_2 <- diag(var_covar_2) # variance of cluster 1  
var_2

sd_2 <- sqrt(var_2) # standard deviation of cluster 1
sd_2




#################
# Principal Components
#################
X.cs=scale(X, scale=T)


# Obtain the PC matrix. 
## First find variance covariance matrix of centered and scaled data
Sxcs=var(X.cs)
Sxcs
####### Looking at the variance covariance matrix, there is a lot of variability in the covariances.
####### However, there are a decent proportion of high covariances (ex. >0.7), so there is high correlation between rows and it is worth it to do PCA

## Then find the eigenvalues (variances of the PC) and the eigenvectors (the axes of the projection)
EP=eigen(Sxcs) 
lambda=EP$values
Proportion.of.variance.of.each.pc=100*(lambda/sum(lambda))
Proportion.of.variance.of.each.pc
cumsum(Proportion.of.variance.of.each.pc)
####### Take note of the cumulative proprtion of variance explained by each PC
####### For my analysis, I'll set the cutoff at 80% i.e. I want my chosen PCs to account for at least 80% of the total variance
####### Thus, I will choose the first 5 PCs, which in total, accounts for 84.7% of the variance  

## Find the eigenvectors
V=EP$vectors
V


# Compute the matrix of principal components scores (the coordinates of the data in the space generated by the eigenvectors )
PC=X.cs%*%V     
head(PC)


# Find the correlation between the components of your PC matrix 
cor(X, PC)
cor(X, PC)[,1:5]   # we are only interested in the first 5 PCs


# biplot
biplot(princomp(X.cs), 
       main = "Biplot of breast cancer data wrt first 2 principal components")


######################
# Perform non-probabilistic k-means clustering using the variables we found significant from the PCA
#####################
datanew2 = data[,c(4,8,9,10,12,24,29)]  # get a dataframe with only the 7 variables that we are interested in
head(datanew2)  # double check

X2 = scale(datanew2, scale=T) # center and scale the data
head(X2)

D2 = dist(X2)   # get the scaled distance matrix
head(D2)



###########
# Non-probabilistic k-means using Lecture 5/5
###########
#### choose initial values  for centers of two clusters 
max(D2)  # tells us the largest distance between rows
D1_2=as.matrix(D2)  #Need this to find rows most distant.

rows.most.distant <- matrix(c(0, 0), ncol = 1) 
# for(i in 1:569){
#   for(j in 1:569){
#     if(D1_2[i, j] == max(D1_2)  ) {
#       rows.most.distant <- matrix(c(i, j), ncol = 1)
#     }
#   }
# }
# rows.most.distant # rows 309 & 79 are the most distant. We set these as initial values. 


# Initial means for the 30 variables in cluster 1 
c1_2=X2[309,] #initial center for cluster 1
c1_2
#Initial means for the 30 variables in cluster 2 
c2_2=X2[79,]  #initial center for cluster 2 
c2_2
#Two index vectors (Z vectors): one  will tell us clusters assigned in last iteration
# Another index vector will tell cluster in the new iteration
# We make the indicator vectors as different as possible to start. 
pastIndicator=569:1   #initial value for z 
indicator=1:569     # past indicator will be compared with new indicator
### note: we initialize this way to get the algorithm started

###### We must iterate until pastIndicator=indicator
## While the two indicator vectors are different, keep going. 
while( sum(pastIndicator!=indicator)!= 0 ) 
{
  pastIndicator=indicator
  
  #distance to current cluster centers
  dc1 = colSums((t(X2)-c1_2)^2)   #distance of each row from mean of cluster 1
  dc2 = colSums((t(X2)-c2_2)^2)   #distance of each row from mean of cluster 2
  dMat = matrix(c(dc1,dc2),,2)
  
  #decide which cluster each point belongs to 
  indicator = max.col(-dMat)
  
  # update the cluster centers
  c1_2 = colMeans(X2[indicator==1,])
  c2_2 = colMeans(X2[indicator==2,])
}

##See the means of each cluster. 
c1_2
c2_2   

## double check that Z_t-1 = Z(t) (convergence)
sum(indicator == pastIndicator) == nrow(datanew) # if we have reached convergence, this should return TRUE

# Use the last Z as indicator for the cluster 
cluster_allocation_2 = indicator
head(cluster_allocation_2)
table(cluster_allocation_2)

#### Get a table with which diagnosis falls in each cluster
#### We can also see if there is any misclassification
table(cluster_allocation_2, diagnosis)  
####### We see from the cluster allocation obtained from using only 7 variables 
####### that We have 162 misclassifications, which means our classification is 71.5% accurate.
####### Percentage wise, it is pretty decent. 
####### So, there is a tradeoff between number of variables included and accuracy of the classification
####### The model you choose depends on how accurate you want your classification to be.  


#### Now we get a table of cluster percentage in each diagnosis
prop_c1_B <- 333/(333+24)
prop_c2_B <- 24/(333+24)
prop_c1_M <- 50/(50+162)
prop_c2_M <- 162/(50+162)
data.frame("cluster" = c(1,2), 'benign'=c(prop_c1_B,prop_c2_B),
           'malignant' = c(prop_c1_M,prop_c2_M))

#### Get a table of diagnosis percentage in each cluster
prop_B_c1 <- 333/(333+50)
prop_M_c1 <- 50/(333+50)
prop_B_c2 <- 24/(24+162)
prop_M_c2 <- 162/(24+162)
data.frame("cluster" = c(1,2), 'benign'=c(prop_B_c1,prop_B_c2),
           'malignant' = c(prop_M_c1,prop_M_c2))

#############
# Newton's method to find the mle of mean radius
#############
# Figuring out which distribution to fit to mean radius
## extract the mean radius from our dataset
data3 = data[,3]
## create a scatterplot to get an idea of the distribution
plot(data3,
     main='Fig IV-1',
     ylab='mean radius')
## get summary statistics for our distribution
summary(data3)
sd(data3)  # standard deviation
var(data3) # variance
## create a potential model to fit our data
plot(rnorm(569, mean=14, sd=3.5), col='red', pch=21,
     main='Fig IV-2',
     ylab='mean radius')


# contour plots to find inital values for mu and sigma^2 
x1=seq(10,20,by=0.1)
x2=seq(10,20,by=0.1)
n=length(data3)
f=matrix(0,nrow=length(x1),ncol=length(x2))
for(i in 1:length(x1)){ 
  for(j in 1:length(x2)) {
    f[i,j]=  -(n/2)*log(x2[j]) -(n/2)*(log(2*pi))-(1/(2*x2[j]))*sum((data3-x1[i])^2)
  }}
contour(x1,x2,f,nlevels=300,xlab="mu",ylab="sigma^2")  
# From the contour plot, we see a critical point around mu=14 and sigma^2=12. So we set those as our inital values.  


#Now numerical optimization
## start defining gradient and hessian
xt=c(100,0)   # this helps us get started 
eps=0.00000001  # tolerance for xtp1-xt 
xtp1=c(14,12)    # vector with initial values for x1(mu) and x2 (sigma^2) 

xHist=matrix(xtp1,2,1) # save history of xt (Stores newton output)
fHist=c()
xHist
### objective function to maximize. The log likelihood of the normal distribution.
### Notice that xtp1[1]= mean mu and xtp1[2] is the variance.
f=-(n/2)*log(xtp1[2])-(n/2)*(log(2*pi))-(1/(2*xtp1[2]))*sum((data3-xtp1[1])^2)

# History of the objective function
fHist=c(f)

while(sum((xtp1-xt)^2)>eps){
  xt=xtp1
  xt
  #compute first and second derivatives
  gradient=as.vector(c((1/xt[2])*(sum(data3-xt[1])) , 
                       -n/(2*xt[2]) + 1/(2*(xt[2]^(2)))*sum((data3 -xt[1])^2)   )) 
  gradient                                
  
  hessian=matrix(c(-n/xt[2],(-1/xt[2]^2)*(sum(data3-xt[1])),(-1/xt[2]^2)*(sum(data3-xt[1])), n/(2*xt[2]^2)-(1/(xt[2]^3))*sum((data3-xt[1])^2) ),ncol=2,nrow=2)
  hessian                                
  ###  
  # compute xtp1 solve(hessian)*gradient=hessian^{-1}*gradient)
  ###   
  xtp1=xt-solve(hessian)%*%gradient # 
  xtp1
  ###   
  #save history
  
  xHist=matrix(c(xHist,xtp1),2)
  xHist
  f=-(n/2)*log(xtp1[2])-(n/2)*(log(2*pi))-(1/(2*xtp1[2]))*sum((data3-xtp1[1])^2)
  
  fHist=c(fHist,f)
}

xHist   #output the estimates of x in each iteration
fHist   #output the f based on estimate of x in each iteration



# Analyze the results
## Get the final gradient and hessian 
gradient                
hessian

# FINAL
fisher = -1*hessian
solve(fisher)

## Get the final estimate
estimate = xHist[,ncol(xHist)]
estimate

## Use the sign of the eigenvalues to find out whether the solution is a minimum, maximum, or saddle point
eigenvals = eigen(hessian)$values 
eigenvals

## Get the variances 
var = diag(solve(-hessian)) 
var

## Then the standard errors 
s.e = sqrt(diag(solve(-hessian)))  
s.e

## Then we can get the 95% CIs for mu and sigma^2  
CI_mu <- estimate[1] + c(-1,1)*1.96*s.e[1] 
CI_mu #confidence interval for mu
CI_sigmasq <- estimate[2] + c(-1,1)*1.96*s.e[2]
CI_sigmasq  #conf int for sigma^2 

#  We are 95% confident that the true parameter mu in the population is between 13.839 and 14.417.  
#  We are 95% confident that the true parameter sigma^2 in the population is between 10.957 and 13.838.   
