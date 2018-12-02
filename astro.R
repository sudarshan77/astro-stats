
library(kedd)
library(mclust)

#xlsx to csv
pulsar_data = read.csv("/J1136+1551_20100106_hi.csv", header=T)
#oulier detection – Using euclidean distance
euclidean_distance = as.matrix( dist( pulsar_data[,-1], upper=T, diag=T) )
n = nrow(euclidean_distance)
outliers <- c()
for(i in 1:nrow(euclidean_distance)){
  pro_obs <- length( which( euclidean_distance[i, -i] <= 0.02 ) ) / n #proportion
  if(pro_obs <= 0.001){
    outliers = union(outliers, i) #appending outlier index
  }
}
# #outliers
length(outliers)

#scatter plot of original data and the data set after removing the outliers located above
plot(pulsar_data[outliers,-1], pch=21, col="red", main="Red Circles indicates Outliers")
points(pulsar_data[-outliers,-1], col='green', pch=21)

# For further analysis Ouliers are removed
pulsar_data = pulsar_data[-outliers,]
#Ploting on_energy - Histogram
hist( pulsar_data[,2], main='on_energy - Histogram', xlab="on_energy", prob=TRUE)
lines( density( pulsar_data[,2] ), col='green')
#Finding the optimal bandwidth for the density estimation for on-pulse energy using
#Epanechnikov kernel using likelihood cross-validation

h.mlcv(x=pulsar_data[,2], lower=0.00001, upper=1, kernel="epanechnikov")

# 0.08916 - optimal bandwith
hist(pulsar_data[,2],prob=TRUE,col="gray",border="white",main=paste("on_pulse Histogram - Epanechnikov kernel of BandWidth - 0.08916"), xlab="on_energy")
lines(density(pulsar_data[,2], bw=0.08916, kernel="epanechnikov"), col='red')
library(mclust)
pulsar_datamodel = Mclust(pulsar_data[,2])
pulsar_datamodel['BIC']



plot(pulsar_datamodel, what="BIC")

plot(pulsar_datamodel,what="density")


pulsar_datamodel$parameters$mean[1]
#off-pulse energies - Mean
mean(pulsar_data[,3])

#Around zero mean is considered to be the null

bootstrap_Clustering = MclustBootstrap(pulsar_datamodel)
bootstrap_summary = summary(bootstrap_Clustering,what='ci')
#confidence interval for the mean of the null component
bootstrap_summary$mean[1:2]
#parameter estimate obtained from pulsar_datamodel
pulsar_datamodel$parameters$mean[1]



#mean of off-pulse energies
mean(pulsar_data[,3])

#off-pulse energies mean is not falling inside confidence interval of mean fornull component
#The mean of off-pulse energies is not same as null component

# Nulling Fraction - Estimation
pulsar_datamodel$parameters$pro[1]

# Nulling Fraction - Confidence Interval
bootstrap_summary$pro[,1]
