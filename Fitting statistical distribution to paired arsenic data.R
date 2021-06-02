rm(list=ls(all=TRUE))

setwd("/Users/yusufjameel/Dropbox/Arsenic_MIT/ICPMS_KIT_analysis")
library(tidyr)
library(MASS)
library(fitdistrplus)
library(gamlss)
library(gamlss.dist)
library(gamlss.add)
library(logspline)
library(ggplot2)
library(geosphere)
library(matrixStats)
library(shape)
library(RColorBrewer)

####################
IC_kit = read.csv("Ahz_ICPMS1.csv", as.is =T)
IC_kit = IC_kit[,c(1:8)]
rem = which(IC_kit$Sample == "1040C")
IC_kit = IC_kit[-rem[1],] ## repeated sample
kit_wide <- spread(IC_kit, kit, lab)
names(kit_wide)[c(7:16)] = paste0("Kit_", "", c(0,10,20,50,100,200,300,500,999,1000))

###########################
kit_0 = subset(kit_wide$Kit_0, kit_wide$Kit_0 >=0)
kit_0[c(which(kit_0 == 0))]= 0.00001

kit_10 = subset(kit_wide$Kit_10, kit_wide$Kit_10 >=0)
kit_10[c(which(kit_10 == 0))]= 0.00001

kit_20 = subset(kit_wide$Kit_20, kit_wide$Kit_20 >=0)
kit_20[c(which(kit_20 == 0))]= 0.00001

kit_50 = subset(kit_wide$Kit_50, kit_wide$Kit_50 >=0)
kit_50[c(which(kit_50 == 0))]= 0.00001

kit_100 = subset(kit_wide$Kit_100, kit_wide$Kit_100 >=0)
kit_100[c(which(kit_100 == 0))]= 0.00001

kit_200 = subset(kit_wide$Kit_200, kit_wide$Kit_200 >=0)
kit_200[c(which(kit_200 == 0))]= 0.00001

kit_300 = subset(kit_wide$Kit_300, kit_wide$Kit_300 >=0)
kit_300[c(which(kit_300 == 0))]= 0.00001

kit_500 = subset(kit_wide$Kit_500, kit_wide$Kit_500 >=0)
kit_500[c(which(kit_500 == 0))]= 0.00001

kit_1000 = subset(kit_wide$Kit_1000, kit_wide$Kit_1000 >=0)
kit_999 = subset(kit_wide$Kit_999, kit_wide$Kit_999 >=0)

kit_1000 = c(kit_1000,kit_999)
kit_1000[c(which(kit_1000 == 0))]= 0.00001


### fitting distributions
  

par(mfrow = c(3,3))

fit_0 <- fitdistr(kit_0, densfun="gamma")  
hist(kit_0, freq = FALSE,  xlim = c(0, quantile(kit_0, 0.9999)))
curve(dgamma(x,  fit_0$estimate[1], fit_0$estimate[2]), from = 0, col = "red", add = TRUE)

fit_10 <- fitdistr(kit_10, densfun="exponential")  
hist(kit_10, freq = FALSE,  xlim = c(0, quantile(kit_10, 0.9999)))
curve(dexp(x,  fit_10$estimate[1]), from = 0, col = "red", add = TRUE)

fit_20 <- fitdistr(kit_20, densfun="gamma")  
hist(kit_20, freq = FALSE,  xlim = c(0, quantile(kit_20, 0.999999)))
curve(dgamma(x,  fit_20$estimate[1], fit_20$estimate[2]), from = 0, col = "red", add = TRUE)

fit_50 <- fitdistr(kit_50, densfun="gamma")  
hist(kit_50, freq = FALSE,  xlim = c(0, quantile(kit_50, 0.99999)))
curve(dgamma(x,  fit_50$estimate[1], fit_50$estimate[2]), from = 0, col = "red", add = TRUE)

fit_100 <- fitdistr(kit_100, densfun="gamma")  
hist(kit_100, freq = FALSE,  xlim = c(0, quantile(kit_100, 0.99999)))
curve(dgamma(x,  fit_100$estimate[1], fit_100$estimate[2]), from = 0, col = "red", add = TRUE)

fit_200 <- fitdistr(kit_200, densfun="gamma")  
hist(kit_200, freq = FALSE,  xlim = c(0, quantile(kit_200, 0.99999)))
curve(dgamma(x,  fit_200$estimate[1], fit_200$estimate[2]), from = 0, col = "red", add = TRUE)

fit_300 <- fitdistr(kit_300, densfun="gamma")  
hist(kit_300, freq = FALSE,  xlim = c(0, quantile(kit_300, 0.99999)))
curve(dgamma(x,  fit_300$estimate[1], fit_300$estimate[2]), from = 0, col = "red", add = TRUE)

fit_500 <- fitdistr(kit_500, densfun="gamma")  
hist(kit_500, freq = FALSE,  xlim = c(0, quantile(kit_500, 0.99999)))
curve(dgamma(x,  fit_500$estimate[1], fit_500$estimate[2]), from = 0, col = "red", add = TRUE)

fit_1000 <- fitdistr(kit_1000, densfun="gamma")  
hist(kit_1000, freq = FALSE,  xlim = c(0, quantile(kit_1000, 0.99999)))
curve(dgamma(x,  fit_1000$estimate[1], fit_1000$estimate[2]), from = 0, col = "red", add = TRUE)

  ### calculate the density of the different kit values associated with the ICPMS measurements
  
  ###once the density is calculated, we normalize the density to calculate the
  ###probabilites of different kit measurements associated with the actual ICP measurements. 
  ###We also calculate the probability of observing red, green or blue color for a given arsenic concnetration



Ahz_ICP = seq(0.5,1000,0.5)

prob_kit0 = matrix(1,length(Ahz_ICP), 1000)
prob_kit10 = matrix(1,length(Ahz_ICP), 1000)
prob_kit20 = matrix(1,length(Ahz_ICP), 1000)
prob_kit50 = matrix(1,length(Ahz_ICP), 1000)
prob_kit100 = matrix(1,length(Ahz_ICP), 1000)
prob_kit200 = matrix(1,length(Ahz_ICP), 1000)
prob_kit300 = matrix(1,length(Ahz_ICP), 1000)
prob_kit500 = matrix(1,length(Ahz_ICP), 1000)
prob_kit1000 = matrix(1,length(Ahz_ICP), 1000)



for(j in 1:1000){
  ICP_random =  runif(length(Ahz_ICP),Ahz_ICP - Ahz_ICP*0.05, Ahz_ICP + Ahz_ICP*0.05 )
  ICP_random[ICP_random <0.5]=0.5
  print(j)
  prob_kit0[,j] = dgamma(ICP_random, fit_0$estimate[1], fit_0$estimate[2])
  prob_kit10[,j] = dexp(ICP_random , fit_10$estimate[1]) 
  prob_kit20[,j] = dgamma(ICP_random , fit_20$estimate[1], fit_20$estimate[2])
  prob_kit50[,j]= dgamma(ICP_random , fit_50$estimate[1], fit_50$estimate[2])
  prob_kit100[,j]= dgamma(ICP_random , fit_100$estimate[1], fit_100$estimate[2])
  prob_kit200[,j] = dgamma(ICP_random , fit_200$estimate[1], fit_200$estimate[2])
  prob_kit300[,j] = dgamma(ICP_random, fit_300$estimate[1], fit_300$estimate[2])
  prob_kit500[,j] = dgamma(ICP_random , fit_500$estimate[1], fit_500$estimate[2])
  prob_kit1000[,j] = dgamma(ICP_random , fit_100$estimate[1], fit_1000$estimate[2])   
}



sum_prob = prob_kit0 + prob_kit10 + prob_kit20 + prob_kit50 + prob_kit100 + prob_kit200 + prob_kit300 + prob_kit500 + prob_kit1000

#############
par(mfrow = c(3,1))

###############plot red green and blue...........
x = matrix(1,length(Ahz_ICP), 100)
y = matrix(1,length(Ahz_ICP), 100)
z = matrix(1,length(Ahz_ICP), 100)

####red placard
for(i in 1:100){
  x[,i] = (prob_kit100[,i]/sum_prob[,i] + prob_kit200[,i]/sum_prob[,i] +
             prob_kit300[,i]/sum_prob[,i] + prob_kit500[,i]/sum_prob[,i]+
             prob_kit1000[,i]/sum_prob[,i])
}

###green placard
y = matrix(1,length(Ahz_ICP), 100)
for(i in 1:100){
  y[,i] = (prob_kit20[,i]/sum_prob[,i] + prob_kit50[,i]/sum_prob[,i])
}

###blue placard
z = matrix(1,length(Ahz_ICP), 100)
for(i in 1:100){
  z[,i] = (prob_kit0[,i]/sum_prob[,i] + prob_kit10[,i]/sum_prob[,i])
}

#####plot the red, green and blue placards
par(mfrow = c(1,1))
matplot(Ahz_ICP, rowMeans(x), type = "l",  axes= F, col=alpha(rgb(1,0,0), 0.9), log = "x",xlab = "Arsenic concentration (ppb)", 
        ylab = "Probability", main = "Probability of color placard assignments", lwd =1)
axis(side = 1, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
axis(side = 2, at = c(0,0.2, 0.4, 0.6,0.8,1), labels = c(0,0.2, 0.4, 0.6,0.8,1))
abline(h = c(0,0.2, 0.4, 0.6,0.8,1), col = '#00000025', lty = "dashed")
abline(v = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
box()
matlines(Ahz_ICP, rowMeans(y), type = "l", col = alpha(rgb(0,1,0), 0.9), lwd =1)
matlines(Ahz_ICP, rowMeans(z), type = "l", col=alpha(rgb(0,0,1), 0.9), lwd=1)


##############plot safe vs. unsafe

x = matrix(1,length(Ahz_ICP), 100)
a = matrix(1,length(Ahz_ICP), 100)

##########unsafe wells
for(i in 1:100){
  x[,i] = (prob_kit100[,i]/sum_prob[,i] + prob_kit200[,i]/sum_prob[,i] +
             prob_kit300[,i]/sum_prob[,i] + prob_kit500[,i]/sum_prob[,i]+
             prob_kit1000[,i]/sum_prob[,i])
}



#############################correct vs. incorrect assignments

matplot(Ahz_ICP[101:2000], 1- rowMeans(x[c(101:2000),]), type = "l",  col=alpha(rgb(1,0,0), 0.9),
        xlim = c(0.5,1000), ylim = c(0,1),axes= F,log = "x",xlab = "Arsenic concentration (ppb)", 
        ylab = "Probability",main = "Probability of incorrect color placard assignments")
axis(side = 1, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
axis(side = 2, at = c(0,0.2, 0.4, 0.6,0.8,1), labels = c(0,0.2, 0.4, 0.6,0.8,1))
abline(h = c(0,0.2, 0.4, 0.6,0.8,1), col = '#00000025', lty = "dashed")
abline(v = c(0.5,10,50,1000), col = '#00000025', lty = "dashed")
box()
matlines(Ahz_ICP[20:100], 1- rowMeans(y[c(20:100),]), type = "l", col = alpha(rgb(0,1,0), 0.9))
matlines(Ahz_ICP[1:20], 1- rowMeans(z[c(1:20),]), type = "l", col = alpha(rgb(0,0,1), 0.9))


######safe wells
a = matrix(1,length(Ahz_ICP), 100)
for(i in 1:100){
  a[,i] = (prob_kit20[,i]/sum_prob[,i] + prob_kit50[,i]/sum_prob[,i]+
             prob_kit0[,i]/sum_prob[,i] + prob_kit10[,i]/sum_prob[,i])
}


matplot(Ahz_ICP, rowMeans(x), type = "l", axes= F, col=alpha(rgb(0,0,0), 0.9), log = "x",xlab = "Arsenic concentration (ppb)", 
        ylab = "Probability",main = "Probability of safe and unsafe assignments")
axis(side = 1, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
axis(side = 2, at = c(0,0.2, 0.4, 0.6,0.8,1), labels = c(0,0.2, 0.4, 0.6,0.8,1))
abline(h = c(0,0.2, 0.4, 0.6,0.8,1), col = '#00000025', lty = "dashed")
abline(v = c(0.5,50,1000), col = '#00000025', lty = "dashed")
box()
matlines(Ahz_ICP, rowMeans(a), type = "l", col = alpha(rgb(0,1,1), 0.9))



