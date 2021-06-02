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

fit_10 <- fitdistr(kit_10, densfun="exponential")  

fit_20 <- fitdistr(kit_20, densfun="gamma")  

fit_50 <- fitdistr(kit_50, densfun="gamma")  

fit_100 <- fitdistr(kit_100, densfun="gamma")  

fit_200 <- fitdistr(kit_200, densfun="gamma")  

fit_300 <- fitdistr(kit_300, densfun="gamma")  

fit_500 <- fitdistr(kit_500, densfun="gamma")  

fit_1000 <- fitdistr(kit_1000, densfun="gamma")  

### calculate the density of the different kit values associated with the ICPMS measurements
###once the density is calculated, we normalize the density to calculate the
###probabilities of different kit measurements associated with the actual ICP measurements. 
###We also calculate the probability of observing red, green or blue color for a given arsenic concentration


Ahz_ICP = read.csv('Araihazar_6000.csv', as.is = T)
rem1 = which(is.na(Ahz_ICP$As_ppb))
rem2 = which(is.na(Ahz_ICP$Latitude))
Ahz_ICP = Ahz_ICP[-c(rem1,rem2),]
prob_kit0 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit10 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit20 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit50 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit100 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit200 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit300 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit500 = matrix(1,nrow(Ahz_ICP), 1000)
prob_kit1000 = matrix(1,nrow(Ahz_ICP), 1000)


############ calculating spatial auto correlation
library(ape)
Ahz_ICP.dists <- as.matrix(dist(cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude)))

Ahz_ICP.inv <- 1/Ahz_ICP.dists
diag(Ahz_ICP.inv) <- 0

Ahz_ICP.inv[which(Ahz_ICP.inv == "Inf")] = 0.001

Ahz_ICP.inv[1:5, 1:5]

Moran.I(Ahz_ICP$As_ppb, Ahz_ICP.inv, na.rm = TRUE)
Moran.I(Ahz_ICP$Random_kit, Ahz_ICP.inv, na.rm = TRUE)

#############################


for(j in 1:1000){
  ICP_random =  rnorm(nrow(Ahz_ICP),Ahz_ICP$As_ppb,Ahz_ICP$As_ppb*0.1)
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



Ahz_ICP$kit_0 = rowMeans(prob_kit0/sum_prob)
Ahz_ICP$kit_10 = rowMeans(prob_kit10/sum_prob)
Ahz_ICP$kit_20  = rowMeans(prob_kit20/sum_prob) 
Ahz_ICP$kit_50  = rowMeans(prob_kit50/sum_prob) 
Ahz_ICP$kit_100  = rowMeans(prob_kit100/sum_prob) 
Ahz_ICP$kit_200  = rowMeans(prob_kit200/sum_prob) 
Ahz_ICP$kit_300  = rowMeans(prob_kit300/sum_prob) 
Ahz_ICP$kit_500  = rowMeans(prob_kit500/sum_prob) 
Ahz_ICP$kit_1000  = rowMeans(prob_kit1000/sum_prob) 



Ahz_ICP$blue = Ahz_ICP$kit_0 + Ahz_ICP$kit_10
Ahz_ICP$green = Ahz_ICP$kit_20 + Ahz_ICP$kit_50
Ahz_ICP$red = Ahz_ICP$kit_100 + Ahz_ICP$kit_200 + Ahz_ICP$kit_300 + Ahz_ICP$kit_500 + Ahz_ICP$kit_1000 

################let's get the kit category for each well using categorical distribution

##draw1
prob_vect = vector()
for ( i in 1: nrow(Ahz_ICP)){
  mulnom_draw = vector() 
  
  prob_vect = c(Ahz_ICP$kit_10[i],Ahz_ICP$kit_0[i],Ahz_ICP$kit_20[i],Ahz_ICP$kit_50[i],Ahz_ICP$kit_100[i],Ahz_ICP$kit_200[i], Ahz_ICP$kit_300[i],Ahz_ICP$kit_500[i],Ahz_ICP$kit_1000[i])
  draw = 1
  mulnom_draw  = rmultinom(1,draw,prob = prob_vect)  ##### get the number of times a kit value would be observed
  
  llk_kit0 = vector()
  llk_kit10= vector()
  llk_kit20= vector()
  llk_kit50= vector()
  llk_kit100= vector()
  llk_kit200= vector()
  llk_kit300= vector()
  llk_kit500= vector()
  llk_kit1000= vector()
  
  
  
  llk_kit0[i] = 0*mulnom_draw[2]
  llk_kit10[i] = 10*mulnom_draw[1]
  llk_kit20[i] = 20*mulnom_draw[3]
  llk_kit50[i] = 50*mulnom_draw[4]
  llk_kit100[i] = 100*mulnom_draw[5]
  llk_kit200[i] = 200*mulnom_draw[6]
  llk_kit300[i] = 300*mulnom_draw[7]
  llk_kit500[i] = 500*mulnom_draw[8]
  llk_kit1000[i] = 1000*mulnom_draw[9]
  
  Ahz_ICP$sim_kit1[i] = (llk_kit0[i]  + llk_kit10[i] + llk_kit20[i] + llk_kit50[i] + llk_kit100[i] + llk_kit300[i] + llk_kit500[i] + llk_kit1000[i] + llk_kit200[i])
}


###draw2
prob_vect = vector()
llk_kit0 = vector()
llk_kit10= vector()
llk_kit20= vector()
llk_kit50= vector()
llk_kit100= vector()
llk_kit200= vector()
llk_kit300= vector()
llk_kit500= vector()
llk_kit1000= vector()
for ( i in 1: nrow(Ahz_ICP)){
  mulnom_draw = vector() 
  
  prob_vect = c(Ahz_ICP$kit_10[i],Ahz_ICP$kit_0[i],Ahz_ICP$kit_20[i],Ahz_ICP$kit_50[i],Ahz_ICP$kit_100[i],Ahz_ICP$kit_200[i], Ahz_ICP$kit_300[i],Ahz_ICP$kit_500[i],Ahz_ICP$kit_1000[i])
  draw = 2
  mulnom_draw  = rmultinom(1,draw,prob = prob_vect)  ##### get the number of times a kit value would be observed
  llk_kit0[i] = 0*mulnom_draw[2]
  llk_kit10[i] = 10*mulnom_draw[1]
  llk_kit20[i] = 20*mulnom_draw[3]
  llk_kit50[i] = 50*mulnom_draw[4]
  llk_kit100[i] = 100*mulnom_draw[5]
  llk_kit200[i] = 200*mulnom_draw[6]
  llk_kit300[i] = 300*mulnom_draw[7]
  llk_kit500[i] = 500*mulnom_draw[8]
  llk_kit1000[i] = 1000*mulnom_draw[9]
  
  Ahz_ICP$sim_kit2[i] = (llk_kit0[i]  + llk_kit10[i] + llk_kit20[i] + llk_kit50[i] + llk_kit100[i] + llk_kit300[i] + llk_kit500[i] + llk_kit1000[i] + llk_kit200[i])
}
##################################IDEAL SWITCHING SCENARIO###################
####################### switching based on actual ICPMS measurements
Ahz_ICP$switch_ICPMS_ideal = Ahz_ICP$As_ppb
Ahz_ICP$Well_ID = seq(1:nrow(Ahz_ICP))
Ahz_ICP$Well_ID_ideal = Ahz_ICP$Well_ID
aaa = Ahz_ICP$Well_ID
rad = 100
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$As_ppb[i]
  if (chk > 0){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    xx =  min(Ahz_ICP$As_ppb[temp]) 
    yy = which.min(Ahz_ICP$As_ppb[temp])
    if (xx < chk){
      Ahz_ICP$switch_ICPMS_ideal[i] = xx 
      aaa[i] = temp[yy]
      print(aaa[i])
    }
  }
}
mean(Ahz_ICP$switch_ICPMS_ideal)
length(which(Ahz_ICP$switch_ICPMS_ideal <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
length(which(Ahz_ICP$switch_ICPMS ==  Ahz_ICP$As_ppb)) ### As exposre did not change

z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}
z
##################################

#################################IDEAL SWITCHING SCENARIO with standard deviation###################
####################### 
ICPMS_SD <- vector()
for (i in 1: nrow(Ahz_ICP)){
  #ICPMS_SD[i] = Ahz_ICP$As_ppb[i] + runif(1,-0.05*Ahz_ICP$As_ppb[i], 0.05*Ahz_ICP$As_ppb[i])
  ICPMS_SD[i] = rnorm(1,Ahz_ICP$As_ppb[i], 0.05*Ahz_ICP$As_ppb[i])
}
ICPMS_SD[ICPMS_SD < 0] = 0
Ahz_ICP$ICPMS_SD = ICPMS_SD
Ahz_ICP$switch_ICPMS_SD = Ahz_ICP$ICPMS_SD
aaa = Ahz_ICP$Well_ID
rad = 100
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$ICPMS_SD[i]
  if (chk > 0){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    xx =  min(Ahz_ICP$ICPMS_SD[temp]) 
    yy = which.min(Ahz_ICP$ICPMS_SD[temp])
    print(xx)
    if (xx < chk){
      Ahz_ICP$switch_ICPMS_SD[i] = xx 
      #aaa[i] = temp[which(Ahz_ICP$ICPMS_SD[temp] == min(Ahz_ICP$ICPMS_SD[temp]))][1]
      aaa[i] = temp[yy]
    }
  }
}
mean(Ahz_ICP$switch_ICPMS_SD)
length(which(Ahz_ICP$switch_ICPMS_SD <  Ahz_ICP$ICPMS_SD))/6595 ### As exposre decreased
length(which(Ahz_ICP$switch_ICPMS_SD ==  Ahz_ICP$ICPMS_SD))/6595 ### As exposre did not change
z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}
z
##################################
####################
####################### optimal switching threshold
optimal_exposure <- vector()
mean_decrease <- vector()


for (a in 1:50){
  print(a)
  Ahz_ICP$switch_ICPMS_optimal = Ahz_ICP$As_ppb
Ahz_ICP$switch_cat2 = 1
Ahz_ICP$switch_cat2[Ahz_ICP$As_ppb > a+10] = 2

kk2 <- Ahz_ICP$As_ppb
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$switch_cat2[i]
  if (chk == 2){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx =  min(Ahz_ICP$switch_cat2[dist_ordered]) 
    yy = which.min(Ahz_ICP$switch_cat2[dist_ordered])
    #print(xx)
    if (xx < chk){
      #Ahz_ICP$switch_cat[i] = xx
      kk2[i] = Ahz_ICP$As_ppb[dist_ordered[yy]]
    }
  }
}
optimal_exposure[a]= mean(kk2)
Ahz_ICP$switch_ICPMS_optimal = kk2
mean_decrease[a] = length(which(Ahz_ICP$switch_ICPMS_optimal <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
}
##
####
mean(Ahz_ICP$switch_ICPMS_optimal)
length(which(Ahz_ICP$switch_ICPMS_optimal <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
length(which(Ahz_ICP$switch_ICPMS_optimal ==  Ahz_ICP$As_ppb))/6595 ### As exposre did not change


##################################
####################### switching based on AAS  categories i.e <10 is blue, <50 is green and >50 is red
Ahz_ICP$switch_cat = 1
Ahz_ICP$switch_cat[Ahz_ICP$As_ppb>10 & Ahz_ICP$As_ppb<50] = 2
Ahz_ICP$switch_cat[Ahz_ICP$As_ppb>50] = 3
aaa = Ahz_ICP$Well_ID
kk2 <- Ahz_ICP$As_ppb
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$switch_cat[i]
  if (chk == 3){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx =  min(Ahz_ICP$switch_cat[dist_ordered]) 
    yy = which.min(Ahz_ICP$switch_cat[dist_ordered])
    print(xx)
    if (xx < chk){
      kk2[i] = Ahz_ICP$As_ppb[dist_ordered[yy]]
      aaa[i] = dist_ordered[yy]
    }
  }
}
mean(kk2)
Ahz_ICP$actual_exp_after_swit_cat = kk2
##
length(which(Ahz_ICP$actual_exp_after_swit_cat <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
length(which(Ahz_ICP$actual_exp_after_swit_cat ==  Ahz_ICP$As_ppb))/6595 ### As exposre did not change



z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}

####################### switching based on AAS  categories i.e  <50 is green and >50 is red
Ahz_ICP$switch_cat2 = 1
Ahz_ICP$switch_cat2[Ahz_ICP$As_ppb > 50] = 2
aaa = Ahz_ICP$Well_ID
kk2 <- Ahz_ICP$As_ppb
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$switch_cat2[i]
  if (chk == 2){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx =  min(Ahz_ICP$switch_cat2[dist_ordered]) 
    yy = which.min(Ahz_ICP$switch_cat2[dist_ordered])
    print(xx)
    if (xx < chk){
      kk2[i] = Ahz_ICP$As_ppb[dist_ordered[yy]]
      aaa[i] = dist_ordered[yy] 
    }
  }
}
mean(kk2)
Ahz_ICP$actual_exp_after_swit_cat2 = kk2
##

####HOW MANY WELLS SWITCHED FROM GREATER THAN 50 LESS THAN 50
length(which(Ahz_ICP$As_ppb >50 & Ahz_ICP$actual_exp_after_swit_cat2<=50))/6595
z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}
z

#####################################################################

##########################################################################
##################################
####################### switching based on simulated kit categories, all category switch except category 0
Ahz_ICP$switch_kit = Ahz_ICP$sim_kit2
aaa = Ahz_ICP$Well_ID
kk1 <- Ahz_ICP$As_ppb
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$sim_kit2[i]
  if (chk > 0){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx =  min(Ahz_ICP$sim_kit2[dist_ordered]) 
    yy = which.min(Ahz_ICP$sim_kit2[dist_ordered])
    print(xx)
    if (xx < chk){
      Ahz_ICP$switch_kit[i] = xx
      kk1[i] = Ahz_ICP$As_ppb[dist_ordered[yy]]
      aaa[i] = dist_ordered[yy] 
    }
  }
}
mean(kk1)
Ahz_ICP$actual_exp_after_swit_kit = kk1

length(which(Ahz_ICP$actual_exp_after_swit_kit <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
length(which(Ahz_ICP$actual_exp_after_swit_kit ==  Ahz_ICP$As_ppb))/6595 ### As exposre did not change
length(which(Ahz_ICP$actual_exp_after_swit_kit > Ahz_ICP$As_ppb))/6595 ### As exposre did not change

z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}
z
#############

####################### switching based on simulated kit categories ABOVE 50 PPB
Ahz_ICP$switch_kit_cat2 = 1
Ahz_ICP$switch_kit_cat2[Ahz_ICP$sim_kit1 > 50] = 2
aaa = Ahz_ICP$Well_ID
kk3 <- Ahz_ICP$As_ppb
zz <- vector()
xx <- vector()
yy <- vector()

for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$switch_kit_cat2[i]
  if (chk == 2){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx[i] =  min(Ahz_ICP$switch_kit_cat2[dist_ordered]) 
    yy[i] = which.min(Ahz_ICP$switch_kit_cat2[dist_ordered])
    print(xx[i])
    if (xx[i] < chk){
      kk3[i] = Ahz_ICP$As_ppb[dist_ordered[yy[i]]]
      aaa[i] = dist_ordered[yy[i]]
    }
  }
}
mean(kk3)
Ahz_ICP$actual_exp_after_swit_kit_cat2 = kk3

length(which(Ahz_ICP$actual_exp_after_swit_kit_cat2 <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
length(which(Ahz_ICP$actual_exp_after_swit_kit_cat2 ==  Ahz_ICP$As_ppb))/6595 ### As exposre did not change
length(which(Ahz_ICP$actual_exp_after_swit_kit_cat2 > Ahz_ICP$As_ppb))/6595 ### As exposre did not change
z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}
z
#############

#######################optimal switching
optimal_exposure_kit <- vector()
mean_decrease_kit <- vector()
mean_increase_kit <- vector()
 for (a in 1:5){
   cat_a = c(10,20,50,100,200)
Ahz_ICP$switch_kit_cat2 = 1
Ahz_ICP$switch_kit_cat2[Ahz_ICP$sim_kit1 > cat_a[a]] = 2
print(a)
kk3 <- Ahz_ICP$As_ppb
zz <- vector()
xx <- vector()
yy <- vector()
aaa = Ahz_ICP$Well_ID
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$switch_kit_cat2[i]

  if (chk == 2){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx[i] =  min(Ahz_ICP$switch_kit_cat[dist_ordered]) 
    yy[i] = which.min(Ahz_ICP$switch_kit_cat[dist_ordered])

    if (xx[i] < chk){
      kk3[i] = Ahz_ICP$As_ppb[dist_ordered[yy[i]]]
      aaa[i] = temp[which(Ahz_ICP$As_ppb[dist_ordered[yy[i]]] == min(Ahz_ICP$As_ppb[dist_ordered[yy[i]]]))][1]
    }
  }
}
optimal_exposure_kit[a] = mean(kk3)
Ahz_ICP$switch_kit_optimal = kk3
mean_decrease_kit[a] = length(which(Ahz_ICP$switch_kit_optimal <  Ahz_ICP$As_ppb))/6595 
mean_increase_kit[a]  = length(which(Ahz_ICP$switch_kit_optimal > Ahz_ICP$As_ppb))/6595
}



length(which(Ahz_ICP$switch_kit_optimal <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
length(which(Ahz_ICP$switch_kit_optimal ==  Ahz_ICP$As_ppb))/6595 ### As exposre did not change
length(which(Ahz_ICP$switch_kit_optimal > Ahz_ICP$As_ppb))/6595 ### As exposre did not change
#############

####################### switching based on simulated kit categories i.e <10 blue, <50 green and <50 red
Ahz_ICP$switch_kit_cat = 1
#Ahz_ICP$switch_kit_cat[Ahz_ICP$sim_kit1 %in% c (20,50)] = 2
Ahz_ICP$switch_kit_cat[Ahz_ICP$sim_kit1 > 50] = 3
aaa = Ahz_ICP$Well_ID
kk3 <- Ahz_ICP$As_ppb
zz <- vector()
xx <- vector()
yy <- vector()

for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$switch_kit_cat[i]
  if (chk == 3){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx[i] =  min(Ahz_ICP$switch_kit_cat[dist_ordered]) 
    yy[i] = which.min(Ahz_ICP$switch_kit_cat[dist_ordered])
    print(xx[i])
    if (xx[i] < chk){
      kk3[i] = Ahz_ICP$As_ppb[dist_ordered[yy[i]]]
      aaa[i] = temp[which(Ahz_ICP$As_ppb[dist_ordered[yy[i]]] == min(Ahz_ICP$As_ppb[dist_ordered[yy[i]]]))][1]
    }
  }
}
mean(kk3)
Ahz_ICP$actual_exp_after_swit_kit_cat = kk3

length(which(Ahz_ICP$actual_exp_after_swit_kit_cat <  Ahz_ICP$As_ppb))/6595 ### As exposre decreased
length(which(Ahz_ICP$actual_exp_after_swit_kit_cat ==  Ahz_ICP$As_ppb))/6595 ### As exposre did not change
length(which(Ahz_ICP$actual_exp_after_swit_kit_cat > Ahz_ICP$As_ppb))/6595 ### As exposre did not change
z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}

#######################################randomise arsenic concentration and check the reduction
############################based upon ICPMS emasurements
Ahz_ICP$Random_As = sample(Ahz_ICP$As_ppb)
aaa = Ahz_ICP$Well_ID
###let's calculate switing basen on this random value

Ahz_ICP$switch_ICPMS_Random = Ahz_ICP$Random_As
rad = 100
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$Random_As[i]
  if (chk > 0){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    xx =  min(Ahz_ICP$Random_As[temp]) 
    yy = which.min(Ahz_ICP$Random_As[temp])
    print(xx)
    if (xx < chk){
      Ahz_ICP$switch_ICPMS_Random[i] = xx 
      aaa[i] = temp[yy]
      
    }
  }
}
mean(Ahz_ICP$switch_ICPMS_Random)

length(which(Ahz_ICP$switch_ICPMS_Random <  Ahz_ICP$Random_As))/6595 ### As exposre decreased
length(which(Ahz_ICP$switch_ICPMS_Random ==  Ahz_ICP$Random_As))/6595 ### As exposre did not change
z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}
z
############################based upon KIT emasurements
df = data.frame(Ahz_ICP$sim_kit2,Ahz_ICP$As_ppb)
df = df[sample(nrow(df), nrow(df)), ]
colnames(df) = c("Random_kit", "Random_As_conc")
###let's calculate switing basen on this random value
aaa = Ahz_ICP$Well_ID
Ahz_ICP$Random_kit = df$Random_kit
Ahz_ICP$Random_As_conc = df$Random_As_conc
kk3 <- Ahz_ICP$Random_As_conc
rad = 100
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$Random_kit[i]
  if (chk > 0){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx[i] =  min(Ahz_ICP$Random_kit[dist_ordered]) 
    yy[i] = which.min(Ahz_ICP$Random_kit[dist_ordered])
    print(xx[i])
    if (xx[i] < chk){
      #Ahz_ICP$switch_kit_cat[i] = xx[i]
      kk3[i] = Ahz_ICP$Random_As_conc[dist_ordered[yy[i]]]
      aaa[i] = dist_ordered[yy[i]]
    }
  }
}
Ahz_ICP$switch_kit_Random = kk3

mean(Ahz_ICP$switch_kit_Random)

length(which(Ahz_ICP$switch_kit_Random <  Ahz_ICP$Random_As_conc))/6595 ### As exposre decreased
length(which(Ahz_ICP$switch_kit_Random ==  Ahz_ICP$Random_As_conc))/6595 ### As exposre did not change
length(which(Ahz_ICP$switch_kit_Random >  Ahz_ICP$Random_As_conc))/6595
z = 0
for (k in 1:nrow(Ahz_ICP)){
  b = aaa[k]
  if(aaa[k] != aaa[b]){
    z = z+1
  }
}
z

############################based upon KIT emasurements - categories

###let's calculate switing basen on this random value
aaa = Ahz_ICP$Well_ID
Ahz_ICP$switch_kit_Random_cat = 1
Ahz_ICP$switch_kit_Random_cat[Ahz_ICP$Random_kit > 50] = 2
rad = 100
kk3 <- Ahz_ICP$Random_As_conc
for ( i in 1:nrow(Ahz_ICP)){
  chk = Ahz_ICP$switch_kit_Random_cat[i]
  if (chk == 2){
    temp =which(distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,] <rad)
    dist_wells = distm(cbind(Ahz_ICP$Longitude[i], Ahz_ICP$Latitude[i]), cbind(Ahz_ICP$Longitude, Ahz_ICP$Latitude))[1,][temp]
    dist_ordered = temp[order(dist_wells)]
    xx[i] =  min(Ahz_ICP$Random_kit[dist_ordered]) 
    yy[i] = which.min(Ahz_ICP$Random_kit[dist_ordered])
    print(xx[i])
    if (xx[i] < chk){
      #Ahz_ICP$switch_kit_cat[i] = xx[i]
      kk3[i] = Ahz_ICP$Random_As_conc[dist_ordered[yy[i]]]
    }
  }
}
Ahz_ICP$switch_kit_Random2 = kk3

mean(Ahz_ICP$switch_kit_Random2)

length(which(Ahz_ICP$switch_kit_Random2 <  Ahz_ICP$Random_As_conc))/6595 ### As exposre decreased
length(which(Ahz_ICP$switch_kit_Random2 ==  Ahz_ICP$Random_As_conc))/6595 ### As exposre did not change
length(which(Ahz_ICP$switch_kit_Random2 >  Ahz_ICP$Random_As_conc))/6595


############### plots 
par(mfrow = c(2,2), mar= c(2,2,1,1))

###ploting of original expousre and after swutichng based upon AAS 
blue = subset(Ahz_ICP, Ahz_ICP$As_ppb <=10)
green = subset(Ahz_ICP, Ahz_ICP$As_ppb>10 & Ahz_ICP$As_ppb<=50)
red = subset(Ahz_ICP, Ahz_ICP$As_ppb>50)
plot( Ahz_ICP$As_ppb,Ahz_ICP$switch_ICPMS_ideal, log ="xy", axes = F, col = "black",xlim = c(1,1000), ylim = c(1,1000),
     xlab="Originial Exposure - ppb",ylab = "Exposure after well switching - ppb",
     main="Spectrometer based ideal switching - Scenario 1")
axis(side = 2, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
abline(h = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(v = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(0,1, lty = "dashed")
box()

###ploting of original expousre and after swutichng based upon color category 
plot(red$As_ppb, red$actual_exp_after_swit_cat, log ="xy", axes = F, col = "red",xlim = c(1,1000), ylim = c(1,1000),
     xlab="Originial Exposure - ppb",ylab = "Exposure after well switching - ppb",
     main="Spectrometer based color placards - Scenario 4")
points(blue$As_ppb, blue$actual_exp_after_swit_cat,  col = "blue")
points(green$As_ppb, green$actual_exp_after_swit_cat,  col = "green")
axis(side = 2, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
abline(h = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(v = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(0,1, lty = "dashed")
box()


###ploting of original expousre and after swutichng based upon kit -kit switching
plot(Ahz_ICP$As_ppb, Ahz_ICP$actual_exp_after_swit_kit, log ="xy", axes = F, col = "black",xlim = c(1,1000), ylim = c(1,1000),
     xlab="Originial Exposure - ppb",ylab = "Exposure after well switching - ppb",
     main="simulated kit categories - Scenario 6")
axis(side = 1, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
axis(side = 2, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
abline(h = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(v = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(0,1, lty = "dashed")
box()

###ploting of original expousre and after swutichng based upon kit -category switching
plot(red$As_ppb, red$actual_exp_after_swit_kit_cat, log ="xy", axes = F, col = "red",xlim = c(1,1000), ylim = c(1,1000),
     xlab="Originial Exposure - ppb",ylab = "Exposure after well switching - ppb",
     main="simulated kit based color categories - Scenario 8")
points(blue$As_ppb, blue$actual_exp_after_swit_kit_cat,  col = "blue")
points(green$As_ppb, green$actual_exp_after_swit_kit_cat,  col = "green")
axis(side = 1, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
axis(side = 2, at = c(0.5,10,50,100,500,1000), labels = c(0,10,50,100,500,1000))
abline(h = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(v = c(0.5,10,50,100,500,1000), col = '#00000025', lty = "dashed")
abline(0,1, lty = "dashed")
box()


###############Plot optimal switching 
par(mfrow = c(1,2))
plot(c(11:60), optimal_exposure,  axes = F, col = "red",xlim = c(10,60), ylim = c(35,65),
     xlab="True arsenic concentration - ppb",ylab = "Mean exposure after well switching - ppb",
     main="Optimal switching concentration - Spectrometer")
axis(side = 1, at = c(10,20,30,40,50,60), labels = c(10,20,30,40,50,60))
axis(side = 2, at = c(35,40,45,50,55,60,65), labels = c(35,40,45,50,55,60,65))
abline(h = c(35,40,45,50,55,60,65), col = '#00000025', lty = "dashed")
abline(v = c(10,20,30,40,50,60), col = '#00000025', lty = "dashed")
box()


plot(c(1:5), optimal_exposure_kit,  axes = F, col = "red",xlim = c(1,5), ylim = c(35,65),
     xlab="Simulated kit categories",ylab = "Mean exposure after well switching - ppb",
     main="Optimal switching concentration - Kit")
axis(side = 1, at = c(1,2,3,4,5), labels = c(2,3,4,5,6))
axis(side = 2, at = c(35,40,45,50,55,60,65), labels = c(35,40,45,50,55,60,65))
abline(h = c(35,40,45,50,55,60,65), col = '#00000025', lty = "dashed")
abline(v = c(1,2,3,4,5), col = '#00000025', lty = "dashed")
box()


