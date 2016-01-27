
#################
### INPUTS
################

#The largest step in precip allowed between time stamps
giantStep<-25

#The largest difference allowed among the wire deltas
deltaThreshold<-0.5

## Low depth range
Depth_Range_Low<--10

##Wetness Threshold
wetnessThreshold<-500


##Number of strain gauges
numStrainGauges<-3


  
  


##Remove any data where one or more strain gauges are still seraching for stability, i.e., coresponding "s" for the strain gauge stability indicating that it is still searching for stability. That is all three strain gagues' stability must read "p" for passed.

#df<-df[(df$stability1=="p" & df$stability2=="p"& df$stability3=="p"),]













################
### UPLOAD DATA
################






##Load raw frequency data

data<-read.table("C:\\Users\\Derek\\Documents\\NEON\\pri_precip\\Golden Data\\priPrecip_gap_goldenData_maxID_13824.csv",header=TRUE,sep=",")



data$time<-as.POSIXct(data$time,format="%Y-%m-%d %H:%M:%OS ",tz="UTC")




###################################################
##CALIBRATION
###################################################
###Convert frequency data into depth data
###################################################

##Gauge 1
c1_1<-0.037648688452734

c2_1<-0.000019848113082

f0_1<-1090.27


##Gauge 2
c1_2<-0.036118604101687

c2_2<-0.000020559340411

f0_2<-1069.67

#Gauge 3
c1_3<-0.041028978386614

c2_3<-0.000017183088270

f0_3<-1055.82



data$DEPTH1<-c1_1*(data$Freq1-f0_1)+(c2_1*((data$Freq1-f0_1)^2))
data$DEPTH2<-c1_2*(data$Freq2-f0_2)+(c2_2*((data$Freq2-f0_2)^2))
data$DEPTH3<-c1_3*(data$Freq3-f0_3)+(c2_3*((data$Freq3-f0_3)^2))

#######################################################################

##create one minute means over the time series
precip.min.means <- aggregate(data[,c(2,3,4,6,7,8)], list(time=cut(data$time,breaks="min")),mean,na.rm = TRUE)

#################################################
##summarize the number of gaps in the data set
data$gapQF<-as.numeric(data$gapQF)
gap.min.means <- aggregate(data$gapQF, list(time=cut(data$time,breaks="min")),sum,na.rm = TRUE)

##If the number of gaps in a one minute average is >=6 then that time stamp is a gap. Thus a vector is added to the precip.min.means 
precip.min.means$gapQF<-ifelse(gap.min.means$x < 6,0,1)

## convert time stamps to posixct
precip.min.means$time<-as.POSIXct(precip.min.means$time)

##Subset the one minute means at five minute increments

#create a vector for the minutes
precip.min.means$min<-format(precip.min.means$time,"%M")


#subset the data set at five minute increments
data<-precip.min.means[precip.min.means$min=="00" |precip.min.means$min=="05"|precip.min.means$min=="10"|precip.min.means$min=="15"|precip.min.means$min=="20"|precip.min.means$min=="25"|precip.min.means$min=="30"|precip.min.means$min=="35"|precip.min.means$min=="40"|precip.min.means$min=="45"|precip.min.means$min=="50"|precip.min.means$min=="55",]






##compute the aveage depth for the strain gauge over the 3 hour period
##First determine the difference between depth measurements. If an NA exists then look two steps back for the difference

#temp vectors for lagged differences
x<-as.numeric(c(NA,diff(data$DEPTH1,lag=1,differences=1)))
y<-as.numeric(c(NA,NA,diff(data$DEPTH1,lag=2,differences=1)))

#final lagged differences for strain gauge 1
data$delta1<-ifelse(!is.na(x),x,y)

#temp vectors for lagged differences
x<-as.numeric(c(NA,diff(data$DEPTH2,lag=1,differences=1)))
y<-as.numeric(c(NA,NA,diff(data$DEPTH2,lag=2,differences=1)))

#final lagged differences for strain gauge 2
data$delta2<-ifelse(!is.na(x),x,y)



#temp vectors for lagged differences
x<-as.numeric(c(NA,diff(data$DEPTH3,lag=1,differences=1)))
y<-as.numeric(c(NA,NA,diff(data$DEPTH3,lag=2,differences=1)))

#final lagged differences for strain gauge 3
data$delta3<-ifelse(!is.na(x),x,y)



##compute the average of all three depths for a given timestamp

del123<-cbind(data[,"delta1"],data[,"delta2"],data[,"delta3"])

data$aveDelta<-apply(del123,MARGIN=1,FUN=mean)








########################
########################
######################
## Precip algorithm that looks at 3 hr chunks of data but process only the last hour
########################
########################
########################

##create indicators for the various flags


precipOut<-c()

dep<-as.numeric(nrow(data))
for(i in 1:((nrow(data)/12)-2)){
  
 

  sindex<-(i*12)-11
  eindex<-24+(i*12)
  x<-data[sindex:eindex,]

  
  ##############################
  ##Wire Weights - Inverse Delta Variance
  ##############################

########################################
#USCRN code chunk - calculateWireweights
########################################

  
  ##Determine numer of complete cases for the delta variance
  d1n<-ifelse(sum(complete.cases(x$delta1))>=sum(complete.cases(x$aveDelta)),sum(complete.cases(x$aveDelta)),sum(complete.cases(x$delta1)))
  
  d2n<-ifelse(sum(complete.cases(x$delta2))>=sum(complete.cases(x$aveDelta)),sum(complete.cases(x$aveDelta)),sum(complete.cases(x$delta2)))
  
  d3n<-ifelse(sum(complete.cases(x$delta3))>=sum(complete.cases(x$aveDelta)),sum(complete.cases(x$aveDelta)),sum(complete.cases(x$delta3)))
  
  
  
  deltaVar1<-1/(sum((x$delta1-x$aveDelta)^2,na.rm=TRUE)/(d1n-1))
  
  deltaVar2<-1/(sum((x$delta2-x$aveDelta)^2,na.rm=TRUE)/(d2n-1))
  
  deltaVar3<-1/(sum((x$delta3-x$aveDelta)^2,na.rm=TRUE)/(d3n-1))
  
  
  ######################################################
  ## Scale the wire weights so that they are out of 1
  ######################################################

####################################
#USCRN code chunk - zeroOutLowWires 
####################################
##Determine if any of the depths from this wire are unreasonably low, i.e., < -10. If TRUE set the delta variance for that wire to zero otherwise keep the previously calcualted delta variance for the wire.
  
##Wire 1
  
if(any(x[,c("DEPTH1")]< Depth_Range_Low,na.rm=TRUE)){
  
  deltaVar1<-0
}else{
  deltaVar1<-deltaVar1}
  
##Wire 2

if(any(x[,c("DEPTH2")]< Depth_Range_Low,na.rm=TRUE)){
  
  deltaVar2<-0
}else{
  deltaVar2<-deltaVar2}


##Wire 3

if(any(x[,c("DEPTH3")]< Depth_Range_Low,na.rm=TRUE)){
  
  deltaVar3<-0
}else{
  deltaVar3<-deltaVar3}

##create indicator based on results
if(any(x[,c("DEPTH1","DEPTH2","DEPTH3")]< Depth_Range_Low,na.rm=TRUE)){
  
lowDepth<-1

}else{
  lowDepth<-0
}




#########################################
#USCRN code chunk - zeroWhereDeltaAnomaly 
#########################################

##Sometimes problems can occur in the first period (i.e., the current hour being processed) where a bad wire returns to normal. If this occurs the wire's weight is set to 0. If the absolute value of a wire's delta is greater than the largest step possible "giantStep" it is defined as a bad wire here and it's weight is set to zero.


##Wire 1

if(any(abs(x[25:36,"delta1"]) > giantStep,na.rm=TRUE)){
  
  deltaVar1<-0
}else{
  deltaVar1<-deltaVar1}

##Wire 2

if(any(abs(x[25:36,"delta2"]) > giantStep,na.rm=TRUE)){
  
  deltaVar2<-0
}else{
  deltaVar2<-deltaVar2}


##Wire 3

if(any(abs(x[25:36,"delta3"]) > giantStep,na.rm=TRUE)){
  
  deltaVar3<-0
}else{
  deltaVar3<-deltaVar3}

##create indicator based on results
if(any(abs(x[25:36,c("delta1","delta2","delta3")]) > giantStep,na.rm=TRUE)){
  
  giant<-1
  
}else{
  giant<-0
}

###################################
#USCRN code chunk - atLeastTwoZeros
###################################

##Check to see if more than two wires are missing/bad if so then set the scaled weight to 0 for each wire so that 0 precip is calculated for the hour. Otherwise calculate the scaled weight for each wire.
if(sum(as.logical(c(deltaVar1==0,deltaVar2==0,deltaVar3==0)))>=2){
  
  ##Scale the weight for wire 1
  
  scaledWireWeight1<-0
  
  ##Scale the weight for wire 2
  
  scaledWireWeight2<-0
  
  ##Scale the weight for wire 3
  
  scaledWireWeight3<-0
  
  ##indicator for missingWireInfoQF
  missingWire<-1
  
  
  }else{

##################################
##USCRN code chunk - normalizeTo1
##################################    

##Determine the total weights
totWeights<-(deltaVar1+deltaVar2+deltaVar3)

##Scale the weight for wire 1

scaledWireWeight1<-deltaVar1/totWeights

##Scale the weight for wire 2

scaledWireWeight2<-deltaVar2/totWeights

##Scale the weight for wire 3

scaledWireWeight3<-deltaVar3/totWeights

##indicator for missingWireInfoQF
missingWire<-0
}

##############################################
#USCRN code chunk - weightedAverageOfDeltas
##############################################

  dep<-cbind(x$time[25:36],x$delta1[25:36],x$delta2[25:36],x$delta3[25:36],(scaledWireWeight1*x$delta1[25:36])+(scaledWireWeight2*x$delta2[25:36])+(scaledWireWeight3*x$delta3[25:36]),rep(lowDepth,times=12),rep(giant,times=12),rep(missingWire,times=12),x$gapQF[25:36])
    
  precipOut<-rbind(precipOut,dep)

  
}



##create column headers
colnames(precipOut)<-c("time","delta1","delta2","delta3","intPrecip","lowDepthQF","exDeltaQF","missingWireInfoQF","gapQF")

##Convert matrix to data frame
precipOut<-as.data.frame(precipOut)

##Convert seconds since epoch timestamp back to posixCT 
precipOut$time<-as.POSIXct(precipOut$time, origin="1970-01-01",tz="UTC")








#########################################
#USCRN code chunk - zeroIfDivergentDeltas
#########################################


##If wire deltas deviate too much from one another set the modified precip value to zero otherwise carry the modified precip value through

zeroIfDivergent(deltas=precipOut[,c("delta1","delta2","delta3")],modifiedPrecip=precipOut$intPrecip,deltaThreshold)



## Convert all NAs in dataframe to zero
precipOut[is.na(precipOut)]<-0


############################################
#USCRN code chunk - zeroIfNonpositiveDelta
#############################################

##Zero out any precip where nonpositive deltas were present from one or more of the strain gauges

zeroIfNonpositiveDelta(deltas=precipOut[,c("delta1","delta2","delta3")],adjustedPrecip=precipOut$adjPrecip)




###############
#Rounding data
###############

##Rounding percipitation values to the tenth decimal place and numbers ending in 5 will be rounded up.

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


##Round the precipitation values
precipOut$roundedPrecip<-round2(precipOut$radjPrecip,n=1)


###############################################
# In place of USCRN code chunk - quantizeValues
###############################################

## Force any value less than 0.2 to 0
precipOut$detectablePrecip<-ifelse(precipOut$roundedPrecip<0.25,0,precipOut$roundedPrecip)


################################################
## NULL if overflowing
################################################

## Force any value less than 0.2 to 0
precipOut$finalPrecip<-ifelse(precipOut$detectablePrecip>1000,0,precipOut$detectablePrecip)
precipOut$overflowQF<-ifelse(precipOut$detectablePrecip>1000,1,0)



plot(precipOut$time,precipOut$finalPrecip)

precipOut$finalcumsum<-cumsum(precipOut$finalPrecip)

##look at the cumulative precip

plot(precipOut$time,precipOut$finalcumsum)

 
library(ggplot2)
foo<-ggplot(precipOut) + geom_point(aes(x = time, y = finalPrecip, colour =  gapQF > 0)) +
  scale_colour_manual(name = 'Gap', values = setNames(c('red','green'),c(T, F))) +
  xlab('Time') + ylab('Precipitation (mm)')

foo

############################################################
############################################################
############################################################
############################################################





#
##
###
####
#####
######
#Zero if Divergent
######
#####
####
###
##
#




## Function to zero out precip values when the deltas among wires differ from one another too much. Threshold is deinfed as deltaThreshold under the inputs. Function finds the minimum and maximum delta for each time stamp and then subtracts the max from the min and assess whether this value is greater than the threshold. If it is then the modified precipitation value is set to zero if not then the modified precipitation value is carried through.

zeroIfDivergent<-function(deltas,modifiedPrecip,deltaThreshold){
  
  maxValue<-apply(deltas,MARGIN=1,FUN=max)
  
  minValue<-apply(deltas,MARGIN=1,FUN=min)
  
  maxMinDiff<-maxValue-minValue
  
  precipOut$adjPrecip<<-ifelse(maxMinDiff>deltaThreshold,0,modifiedPrecip)
  
  precipOut$gaugeNoiseQF<<-ifelse(maxMinDiff>deltaThreshold,1,0)
}




#
##
###
####
#####
######
#Zero precip if nonpositive delta
######
#####
####
###
##
#



## Function that assess whether any of the delta's from the three wires is not positive. If so then the adjPrecip value is set to zero if not then the adjPrecip value is carried through.


zeroIfNonpositiveDelta<-function(deltas,adjustedPrecip){
  
  precip<-c()
  results<-c()
  wireNoise<-c()
  
  ##Convert any NAs to 0
  deltas[is.na(deltas)]<-0
  
  for(u in 1:nrow(precipOut)){
    
    ##If any of the deltas are non-positive then set the precip to 0
    if(any(deltas[u,]<0)){
      
      precip[u]<-0
      
      wireNoise[u]<-1
      
    }
    
    ##Otherwise if the deltas are all >=0 then carry through the precip value
    else{
      
      precip[u]<-adjustedPrecip[u]
      
      wireNoise[u]<-0
    }
    
    results<-c(precip)
    wireSum<-c(wireNoise)
  
  }
  
  ##write the results to the precipOut data frame
  precipOut$radjPrecip<<-results
  precipOut$wireNoiseQF<<-wireSum
}













































##add in WAVGPRCP
dataSub<-data[,c("time","WAVGPRCP","PAIRWISEPRCP","DEPTH1","DEPTH2","DEPTH3")]

finalPrecipCheck<-merge(precipOut,dataSub,by="time")


write.table(finalPrecipCheck,"C:\\Users\\dsmith\\Documents\\FIU\\Precipitation\\Primary Precip Gauge\\USCRN data\\NEON_and_USCRN_Precip_BG_KY_May_2013.csv",sep=",")


head(precipOut)
newfin<-precipOut[,c(1:5,19,21,22,23)]
head(newfin)

colnames(newfin)<-c("time","wetness","delta1","delta2","delta3","NEON_calc_no_wetness","NEON_calc_wetness","WAVGPRCP","PAIRWISEPRCP")
###################
######################
#######################
#########################
###########################
############################
##load results for ggplot

precipOut<-read.table("C:\\Users\\dsmith\\Documents\\FIU\\Precipitation\\Primary Precip Gauge\\USCRN data\\finalDepth_USCRN_and_neon_last_hr.csv",sep=",",header=TRUE)

precipOut$time<-as.POSIXct(precipOut$time,format="%Y-%m-%d %H:%M:%S",tz="UTC")


##Subsetting precip data to be used with GGplot
#options "time","finalPrecipReg","finalPrecipRegNonNeg","finalPrecipWet","finalPrecipWetNonNeg","WAVGPRCP"
#

subsetOfPrecipOut<-precipOut[,c("time","finalPrecipReg","finalPrecipRegNonNeg")]

##Convert NA to zero
precipOut[is.na(precipOut)]<-0


library(reshape2)
library(ggplot2)
library(RColorBrewer)
##Melt data for use in ggplot

meltedFinalPrecipCheck<-melt(precipOut,id.vars = c("time"))


##options to plot
# "WAVGPRCP","finalPrecipReg","finalPrecipRegNonNeg","finalPrecipWet","finalPrecipWetNonNeg"

regNonNeg<-c("finalPrecipRegNonNeg")
reg<-c("finalPrecipReg")
wet<-c("finalPrecipWet")
wetNonNeg<-c("finalPrecipWetNonNeg")
uscrnWavgPrcip<-c("WAVGPRCP")



###Plot in ggplot
figall<-ggplot()+
  geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==uscrnWavgPrcip,],aes(x=time,y=value,group=variable,colour=variable))+
  geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wetNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==reg,],aes(x=time,y=value,group=variable,colour=variable))+
 

  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==regNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wet,],aes(x=time,y=value,group=variable,colour=variable))+
  #scale_x_datetime(breaks = "5 hour",labels = date_format("%H"), minor_breaks="day")+
 scale_colour_brewer(name="Precipitation Algorithm ",palette="Set1",labels=c("NEON calc (using wetness sensor data)", "NEON calc (w/out wetness sesnor data)","USCRN WAVGPRCP"))+
                        #,"2 (Mean)","2 (RMSE)","5 (Mean)","5 (RMSE)"))+
  labs(x="Time (UTC)", y="Precipitation (mm)")+
  
  theme(title = element_text(size = 14, colour = 'black'))+
  theme_bw()

figall



figwavg<-ggplot()+
  geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==uscrnWavgPrcip,],aes(x=time,y=value), col="green3")+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wetNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==reg,],aes(x=time,y=value,group=variable,colour=variable))+
  
  
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==regNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wet,],aes(x=time,y=value,group=variable,colour=variable))+
  #scale_x_datetime(breaks = "5 hour",labels = date_format("%H"), minor_breaks="day")+
 # scale_colour_brewer(name="Precipitation Algorithm ",palette="Set1",labels=c("NEON calc (using wetness sensor data)", "NEON calc (w/out wetness sesnor data)","USCRN WAVGPRCP"))+
  #,"2 (Mean)","2 (RMSE)","5 (Mean)","5 (RMSE)"))+
  labs(x="Time (UTC)", y="Precipitation (mm)")+
  
  theme(title = element_text(size = 14, colour = 'black'))+
  theme_bw()

figwavg




figwet<-ggplot()+
  geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wetNonNeg,],aes(x=time,y=value), col="firebrick3")+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wetNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==reg,],aes(x=time,y=value,group=variable,colour=variable))+
  
  
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==regNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wet,],aes(x=time,y=value,group=variable,colour=variable))+
  #scale_x_datetime(breaks = "5 hour",labels = date_format("%H"), minor_breaks="day")+
  # scale_colour_brewer(name="Precipitation Algorithm ",palette="Set1",labels=c("NEON calc (using wetness sensor data)", "NEON calc (w/out wetness sesnor data)","USCRN WAVGPRCP"))+
  #,"2 (Mean)","2 (RMSE)","5 (Mean)","5 (RMSE)"))+
  labs(x="Time (UTC)", y="Precipitation (mm)")+
  
  theme(title = element_text(size = 14, colour = 'black'))+
  theme_bw()

figwet




figReg<-ggplot()+
  geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==reg,],aes(x=time,y=value), col="dodgerblue3")+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wetNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==reg,],aes(x=time,y=value,group=variable,colour=variable))+
  
  
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==regNonNeg,],aes(x=time,y=value,group=variable,colour=variable))+
  #geom_point(data = meltedFinalPrecipCheck[meltedFinalPrecipCheck$variable==wet,],aes(x=time,y=value,group=variable,colour=variable))+
  #scale_x_datetime(breaks = "5 hour",labels = date_format("%H"), minor_breaks="day")+
  # scale_colour_brewer(name="Precipitation Algorithm ",palette="Set1",labels=c("NEON calc (using wetness sensor data)", "NEON calc (w/out wetness sesnor data)","USCRN WAVGPRCP"))+
  #,"2 (Mean)","2 (RMSE)","5 (Mean)","5 (RMSE)"))+
  labs(x="Time (UTC)", y="Precipitation (mm)")+
  
  theme(title = element_text(size = 14, colour = 'black'))+
  theme_bw()

figReg

multiplot(figall,figwavg,figwet,figReg)


summary(finalPrecipCheck)
plot(finalPrecipCheck$time,finalPrecipCheck$finalPrecipReg)
points(finalPrecipCheck$time,finalPrecipCheck$finalPrecipRegNonNeg,col="yellow")

summary(precipOut)

precipOut[precipOut==0]<-NA
stats(precipOut)

library(fields)
#check which rows are NA
which(is.na(finalPrecipCheck), arr.ind=TRUE)
finalPrecipCheck[6840:6855,]
##########
##PLOTS
##########

##Plot WAVGPRCP
plot(data$time,data$WAVGPRCP,col="black")

##Add points from the USCRN pairwise algoirthm
points(data$time,data$PAIRWISEPRCP,col="red")

##Add points from calculating the precp for the last hour only and using three hours to compute the weights
points(precipOut$time,precipOut$precipOut$finalPrecipReg,col="purple")


##Add points from calculating the precp for three hour chunks
points(calcDepth_3hr$time,calcDepth_3hr$precip,col="orange")

##Add points from calculating the precp by using three hours of data but only processing the last point and then increment the window by one
points(calcDepth_inc$time,calcDepth_inc$precip,col="green")











#
##
###
####
#####
######
#Remove Giant Steps
######
#####
####
###
##
#




###Function to remove giant steps in the calculated precipitation values. Compares calcuated precipitation values to the largest step that is feasible over the averaging period.Inputs are,
# - giantStep threshold value
# - precip = the column of of precipOut that contains the initial calculated precipitation values

removeGiantSteps<-function(precip,giantStep){

precipOut$sanePrecip<<-ifelse(abs(precip)>giantStep,0,precip)

}




#
##
###
####
#####
######
#Zero if Divergent
######
#####
####
###
##
#




## Function to zero out precip values when the deltas among wires differ from one another too much. Threshold is deinfed as deltaThreshold under the inputs. Function finds the minimum and maximum delta for each time stamp and then subtracts the max from the min and assess whether this value is greater than the threshold. If it is then the modified precipitation value is set to zero if not then the modified precipitation value is carried through.

zeroIfDivergent<-function(deltas,modifiedPrecip,deltaThreshold){
  
  maxValue<-apply(deltas,MARGIN=1,FUN=max)
  
  minValue<-apply(deltas,MARGIN=1,FUN=min)
  
  maxMinDiff<-maxValue-minValue
  
  precipOut$adjPrecip<<-ifelse(maxMinDiff>deltaThreshold,0,modifiedPrecip)
}




#
##
###
####
#####
######
#Zero precip if nonpositive delta
######
#####
####
###
##
#



## Function that assess whether any of the delta's from the three wires is not positive. If so then the adjPrecip value is set to zero if not then the adjPrecip value is carried through.


zeroIfNonpositiveDelta<-function(deltas,adjustedPrecip){
  
  precip<-c()
  results<-c()
  
  ##Convert any NAs to 0
  deltas[is.na(deltas)]<-0
  
  for(u in 1:nrow(precipOut)){
    
    ##If any of the deltas are non-positive then set the precip to 0
    if(any(deltas[u,]<0)){
      
      precip[u]<-0
      
    }
    
    ##Otherwise if the deltas are all >=0 then carry through the precip value
    else{
      
      precip[u]<-adjustedPrecip[u]
      
    }
    
    results<-c(precip)
  }
  
  ##write the results to the precipOut data frame
  precipOut$radjPrecip<<-results
}




#########
##############
###################
######################
##Wetness Sensor Data
#######################
##################
################
###########
##########
######

res<-c()
outpre<-c()

wetnessSensor<-function(calcPrecip,wetness){
  for(n in 1:length(calcPrecip)){
    
    if(wetness[n]<wetnessThreshold){
      
      outpre[n]<-calcPrecip[n]
    }else{
      outpre[n]<-0
    }
    res<-c(outpre)
  }
  
  precipOut$wetAdjPrecip<<-res
  
}



##Check wetness sensor data
wetnessSensor(calcPrecip=precipOut$modPrecip,wetness=precipOut$wetness)





#
##
###
####
#####
######
#Get nonnegative precip (This chunk of code from USCRN doesn't make sense)
######
#####
####
###
##
#
head(precipOut)
precipOut$wetAdjPrecip
precipOut$modPrecip
##This chunk adjusts precip values so those that removed for having negative deltas do no affect the calculation of the total precip.
# 
 getNonnegativePrecip()




getNonnegativePrecip<-function(){
  
  totPos<-sum(precipOut$wetAdjPrecip[precipOut$wetAdjPrecip >0])
  
  totAll<-sum(precipOut$wetAdjPrecip )
  
  corrfact<-totAll/totPos
  
  finPrecip<-precipOut$wetAdjPrecip *corrfact
  
  precipOut$nonnegativePrecip_Wetness<<-finPrecip
}

################################### 

########################################################
##########################################################
###NEED to make this function have a 3 hour sliding window
##########################################################
##########################################################

  
  


  pre<-precipOut$wetAdjPrecip
  
  getNonnegativePrecip<-function(pre){
    
    negCorPrecipHr<-as.numeric(length(pre))
    negCorPrecipOut<-rep(NA,times=24)
    
    for(a in 1:((length(pre)/12)-2)){
      
      sindex<-(a*12)-11
      eindex<-24+(a*12)
      x<-pre[sindex:eindex]
      
      if(sum(x[x>0])==0){
        
        negCorPrecipHr<-rep(0,times=12)
        
      }else{
      
      totPos<-sum(x[x>0])
      
      totAll<-sum(x)
      
      corFact<-totAll/totPos
      
      negCorPrecip<-x *corFact
      
      negCorPrecipHr<-negCorPrecip[25:36]
      
      
      }
      
      negCorPrecipOut<-c(negCorPrecipOut,negCorPrecipHr)
    }
    
    negCorPrecipOut<-ifelse(negCorPrecipOut<0,0,negCorPrecipOut)
    

    precipOut$nonnegativePrecip_with_out_wetness_sensor<<-negCorPrecipOut
  }
    
  
  getNonnegativePrecip(pre=precipOut$radjPrecip)
  
    length(precipOut$radjPrecip)
summary(precipOut)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  



##Rounding percipitation values to the tenth decimal place and numbers ending in 5 will be rounded up

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


##Round the precipitation values
precipOut$roundedPrecip<-round2(precipOut$radjPrecip,n=1)


## Force any value less than 0.2 to 0
precipOut$finalPrecip<-ifelse(precipOut$roundedPrecip<0.2,0,precipOut$roundedPrecip)





#####################
#####################
### END
#####################
#####################





lasthrprecip<-precipOut


lasthr<-lasthrprecip[,c("time","finalPrecip")]
incby1<-precipOut[,c("time","finalPrecip")]

final<-merge(lasthr,incby1,by="time",all=TRUE)
colnames(final)<-c("time","lasthr","incby1")
points(lasthrprecip$time,lasthrprecip$finalPrecip)
plot(precipOut$time,precipOut$finalPrecip,col="red")
plot(final$time,final$lasthr)
points(final$time,final$incby1,col="red")

tail(lasthrprecip,n=500)
##########
##PLOTS
##########

##Plot WAVGPRCP
plot(data$time,data$WAVGPRCP,col="black")

##Add points from the USCRN pairwise algoirthm
points(data$time,data$PAIRWISEPRCP,col="red")

##Add points from calculating the precp for the last hour only and using three hours to compute the weights
points(calcDepth_1hr$time,calcDepth_1hr$precip,col="purple")


##Add points from calculating the precp for three hour chunks
points(calcDepth_3hr$time,calcDepth_3hr$precip,col="orange")

##Add points from calculating the precp by using three hours of data but only processing the last point and then increment the window by one
points(calcDepth_inc$time,calcDepth_inc$precip,col="green")


##Final calculated precipitation value for 3 hr of data but onl last hour processed
points(precipOut$time,precipOut$finalPrecip, col="red")

##Final calculated precipitation value for 3 hr of data but onl last hour processed
points(precipOut$time,precipOut$roundedPrecip, col="blue")

abline(h=0.2)
precipOut$radjPrecip[5000:5500]



inc1precip<-precipOut
lasthrprecip<-precipOut
last3hrprecip<-precipOut


##Plot WAVGPRCP
plot(data$time,data$WAVGPRCP,col="black")

##Final calculated precipitation value for 3 hr of data but inc every timestamp
plot(inc1precip$time,inc1precip$finalPrecip, col="green")

##Final calculated precipitation value for 3 hr of data
points(precipOut$time,precipOut$finalPrecip, col="red")

##Final calculated precipitation value for 3 hr of data
points(precipOut$time,precipOut$finalPrecipWet, col="blue")

##Final calculated precipitation value for 3 hr of data but onl last hour processed
points(lasthrprecip$time,lasthrprecip$finalPrecip, col="red")


##Round the precipitation values
last3hrprecip$roundedPrecip<-round2(last3hrprecip$radjPrecip,n=1)


## Force any value less than 0.2 to 0
last3hrprecip$finalPrecip<-ifelse(last3hrprecip$roundedPrecip<0.2,0,last3hrprecip$roundedPrecip)

hist(precipOut$finalPrecip[!precipOut$finalPrecip==0],breaks=len)
hist(data$WAVGPRCP[!data$WAVGPRCP==0],breaks=len)

len<-seq(from=0, to=8,by=0.1)
head(precipOut$finalPrecip[!precipOut$finalPrecip==0],n=100)

length(data$WAVGPRCP[!data$WAVGPRCP==0])


length(precipOut$finalPrecip[!precipOut$finalPrecip==0])


