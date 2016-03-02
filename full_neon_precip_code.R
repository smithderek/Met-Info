###############################################################################################################
## Inputs:
################################################################################################################
#location of data files
loc<-"C:\\Users\\dsmith\\Documents\\FIU\\Git\\Met-Info\\pri_precip\\Golden Data\\TarnishedData_2016_02_10"
################################################################
### Thresholds:
#The largest step in precip allowed between time stamps
giantStep<-25

#The largest difference allowed among the wire deltas
deltaThreshold<-0.5

## Low depth range
Depth_Range_Low<--10

################################################################
###Calibration info:
# Using calibration coefficients from maximo ID #13824

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
###############################################################################################################
###############################################################################################################


#################################################################
## Load files and organise into correct format for the algorithm
##################################################################



##list paths to the functions
pathnames <- list.files(pattern="[.]csv$", path=loc, full.names=TRUE)

##source the functions from the given path
dataAsList<-lapply(pathnames, FUN=read.csv,stringsAsFactors=FALSE,header=TRUE)



multmerge <- function(mypath){
  filenames=list.files(pattern="[.]csv$", path=mypath, full.names=TRUE)
  datalist = lapply(filenames, function(x){read.csv(file=x,header=TRUE,stringsAsFactors=FALSE)})
  Reduce(function(x,y) {merge(x,y,by="time")}, datalist)
}



data<-multmerge(mypath=loc)


#Check if headers are what is expected in the algorithm
if(all(colnames(data)==c("time" , "Freq1","stab1","Freq2","stab2","Freq3","stab3","orificeHeater","inletTemp","internTemp"))==FALSE)
{print("column names are inconsistent")}


#convert back from CI timestamp to R timestamp
data$time<-as.POSIXct(data$time,format="%d-%b-%y %I.%M.%OS %p",tz="UTC")

#in merging all of the data the time stamps are out of order
#sort the merged data sets by time since they are
data<-data[order(data$time),]


######################################################################################################
## Calibration
######################################################################################################

###########################################
## Convert to frequency values
#########################################

data$DEPTH1<-c1_1*(data$Freq1-f0_1)+(c2_1*((data$Freq1-f0_1)^2))
data$DEPTH2<-c1_2*(data$Freq2-f0_2)+(c2_2*((data$Freq2-f0_2)^2))
data$DEPTH3<-c1_3*(data$Freq3-f0_3)+(c2_3*((data$Freq3-f0_3)^2))

######################################################################################################
## Calibration end
######################################################################################################








#######################################################################################################
### Precip Calc start
#######################################################################################################
##########################################################
##########################################################


##Check stability info 
data$Freq1<-ifelse(data$stab1=="1",data$Freq1,NA)
data$Freq1<-ifelse(data$stab1=="1",data$Freq1,NA)
data$Freq1<-ifelse(data$stab1=="1",data$Freq1,NA)

##Fail
data$wireStafailQF1<-ifelse(data$stab1=="1" | is.na(data$stab1) ,0,100)
data$wireStafailQF2<-ifelse(data$stab2=="1"| is.na(data$stab2),0,100)
data$wireStafailQF3<-ifelse(data$stab3=="1"| is.na(data$stab3),0,100)

##Pass
data$wireStaPassQF1<-ifelse(data$stab1=="0"| is.na(data$stab1),0,100)
data$wireStaPassQF2<-ifelse(data$stab2=="0"| is.na(data$stab2),0,100)
data$wireStaPassQF3<-ifelse(data$stab3=="0"| is.na(data$stab3),0,100)

##Missing
data$wireStaNAQF1<-ifelse(is.na(data$stab1),100,0)
data$wireStaNAQF2<-ifelse(is.na(data$stab2),100,0)
data$wireStaNAQF3<-ifelse(is.na(data$stab3),100,0)


## Need to average the QM over the one minute period
##############################################################################
##############################################################################
## Orifice heater info
##############################################################################


data$inletHeat1<-ifelse(data$orificeHeater=="100",100,0)

data$inletHeat2<-ifelse(data$orificeHeater=="110",100,0)

data$inletHeat3<-ifelse(data$orificeHeater=="111",100,0)

data$inletHeatNA<-ifelse(is.na(data$orificeHeater),100,0)



inletHeatQMs <- aggregate(data[,c("inletHeat1","inletHeat2","inletHeat3","inletHeatNA")], list(time=cut(data$time,breaks="5 min")),FUN="mean",na.rm = TRUE)

inletHeatQMs$time<-as.POSIXct(inletHeatQMs$time,tz="UTC")

#################################################################################################
#################################################################################################

###################################################################################################
## Aggregate the rpecip data into one minute averages that will be captured every five minutes
####################################################################################################


##create one minute means over the time series
precip.min.means <- aggregate(data[,c("Freq1","Freq2","Freq3","DEPTH1","DEPTH2","DEPTH3","wireStafailQF1","wireStafailQF2","wireStafailQF3","wireStaPassQF1","wireStaPassQF2","wireStaPassQF3","wireStaNAQF1","wireStaNAQF2","wireStaNAQF3")], list(time=cut(data$time,breaks="min")),mean,na.rm = TRUE)

#write.table(data,"C:\\Users\\dsmith\\Desktop\\data_test.csv",sep=",",row.names=FALSE)

#################################################
## GAPS
#################################################
## create indicator for gaps in the data set
#################################################
#
#indicator for a gap
gap<-rep(1,nrow(data))

#indicator for no gap
noGap<-rep(0,nrow(data))

#logical test for gaps
gapTest<-as.logical(is.na(data$Freq1),is.na(data$Freq2),is.na(data$Freq3))

#set true gaps indicator to 1 and 0 otherwise
data$gapQF<-ifelse(gapTest==TRUE,gap,noGap)


##summarize the number of gaps in the data set
data$gapQF<-as.numeric(data$gapQF)
gap.min.means <- aggregate(data$gapQF, list(time=cut(data$time,breaks="min")),sum,na.rm = TRUE)

##If the number of gaps in a one minute average is >=6 then that time stamp is a gap. Thus a vector is added to the precip.min.means 
precip.min.means$gapQF<-ifelse(gap.min.means$x < 6,0,1)

## convert time stamps to posixct
precip.min.means$time<-as.POSIXct(precip.min.means$time,tz="UTC")

##Subset the one minute means at five minute increments

#create a vector for the minutes
precip.min.means$min<-format(precip.min.means$time,"%M")


#subset the data set at five minute increments
data<-precip.min.means[precip.min.means$min=="00" |precip.min.means$min=="05"|precip.min.means$min=="10"|precip.min.means$min=="15"|precip.min.means$min=="20"|precip.min.means$min=="25"|precip.min.means$min=="30"|precip.min.means$min=="35"|precip.min.means$min=="40"|precip.min.means$min=="45"|precip.min.means$min=="50"|precip.min.means$min=="55",]



##add in the 5 minute averages of the orifice heater flag
data<-merge(data,inletHeatQMs,by="time")



##compute the aveage depth for the strain gauge over the 3 hour period
##First determine the difference between depth measurements. If an NA exists then look two steps back for the difference if that doesn't exsit then precip is set to NA and if measurements are missing from two or more gauges the priorDepthQF is set to 1.

#temp vectors for lagged differences
x<-as.numeric(c(NA,diff(data$DEPTH1,lag=1,differences=1)))
y<-as.numeric(c(NA,NA,diff(data$DEPTH1,lag=2,differences=1)))
z<-data$DEPTH1

#final lagged differences for strain gauge 1
data$delta1<-ifelse(!is.na(x),x,y)

data$priorDepth1<-ifelse(is.na(z),0,ifelse(!is.na(x),0,ifelse(!is.na(y),0,1)))

#temp vectors for lagged differences
x<-as.numeric(c(NA,diff(data$DEPTH2,lag=1,differences=1)))
y<-as.numeric(c(NA,NA,diff(data$DEPTH2,lag=2,differences=1)))
z<-data$DEPTH2

#final lagged differences for strain gauge 2
data$delta2<-ifelse(!is.na(x),x,y)

data$priorDepth2<-ifelse(is.na(z),0,ifelse(!is.na(x),0,ifelse(!is.na(y),0,1)))

#temp vectors for lagged differences
x<-as.numeric(c(NA,diff(data$DEPTH3,lag=1,differences=1)))
y<-as.numeric(c(NA,NA,diff(data$DEPTH3,lag=2,differences=1)))
z<-data$DEPTH3

#final lagged differences for strain gauge 3
data$delta3<-ifelse(!is.na(x),x,y)

data$priorDepth3<-ifelse(is.na(z),0,ifelse(!is.na(x),0,ifelse(!is.na(y),0,1)))

##create the priorDepthQF that is set to one if data is missing for 2 or more of the strain gauges

#create a vector with the sums of the priorDepth results
sumOfpriorDepths<-apply(data[,c("priorDepth1","priorDepth2", "priorDepth3")],MARGIN=1,FUN=sum)

#create the priorDepthQF
data$priorDepthQF<-ifelse(sumOfpriorDepths >= 2,1,0)



##compute the average of all three depths for a given timestamp

del123<-cbind(data[,"delta1"],data[,"delta2"],data[,"delta3"])

data$aveDelta<-apply(del123,MARGIN=1,FUN=mean)




################################################
## NULL if overflowing
################################################

## 

depthInfo<-apply(data[,c("DEPTH1","DEPTH2", "DEPTH3")],MARGIN=1,FUN=sum)

#check to see if any of the depth averages indicate the gauge was overflowing
depth1OF<-ifelse(data$DEPTH1>=1000,1,0)
depth2OF<-ifelse(data$DEPTH2>=1000,1,0)
depth3OF<-ifelse(data$DEPTH3>=1000,1,0)

##combine the results
depthOF<-cbind(depth1OF,depth2OF,depth3OF)
depthOFInfo<-apply(depthOF,MARGIN=1,FUN=sum)
##set the precip to 0 if the gauge was overflowing and set the overflowQF to 1
data$overflowQF<-ifelse(depthOFInfo>=1,1,0)















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
  
  dep<-cbind(x$time[25:36],x$delta1[25:36],x$delta2[25:36],x$delta3[25:36],(scaledWireWeight1*x$delta1[25:36])+(scaledWireWeight2*x$delta2[25:36])+(scaledWireWeight3*x$delta3[25:36]),rep(lowDepth,times=12),rep(giant,times=12),rep(missingWire,times=12),x$gapQF[25:36],x$priorDepthQF[25:36],x$overflowQF[25:36],x$inletHeat1[25:36],x$inletHeat2[25:36],x$inletHeat3[25:36],x$inletHeatNA[25:36],x$wireStafailQF1[25:36],x$wireStafailQF2[25:36],x$wireStafailQF3[25:36],x$wireStaPassQF1[25:36],x$wireStaPassQF2[25:36],x$wireStaPassQF3[25:36],x$wireStaNAQF1[25:36],x$wireStaNAQF2[25:36],x$wireStaNAQF3[25:36])
  
  precipOut<-rbind(precipOut,dep)
  
  
}



##create column headers
colnames(precipOut)<-c("time","delta1","delta2","delta3","intPrecip","lowDepthQF","exDeltaQF","missingWireInfoQF","gapQF", "priorDepthQF","overflowQF","inletHeaters1QM","inletHeaters2QM","inletHeaters3QM","inletHeatersNAQM","wire1StabilityFailQM","wire2StabilityFailQM","wire3StabilityFailQM","wire1StabilityPassQM","wire2StabilityPassQM","wire3StabilityPassQM","wire1StabilityNAQM","wire2StabilityNAQM","wire3StabilityNAQM")

##Convert matrix to data frame
precipOut<-as.data.frame(precipOut)


##Convert seconds since epoch timestamp back to posixCT 
precipOut$time<-as.POSIXct(precipOut$time, origin="1970-01-01",tz="UTC")

## Convert all NAs in dataframe to zero
precipOut[is.na(precipOut)]<-0





#########################################
#USCRN code chunk - zeroIfDivergentDeltas
#########################################


## Function to zero out precip values when the deltas among wires differ from one another too much. Threshold is deinfed as deltaThreshold under the inputs. Function finds the minimum and maximum delta for each time stamp and then subtracts the max from the min and assess whether this value is greater than the threshold. If it is then the modified precipitation value is set to zero if not then the modified precipitation value is carried through.

zeroIfDivergent<-function(deltas,modifiedPrecip,deltaThreshold){
  
  maxValue<-apply(deltas,MARGIN=1,FUN=max)
  
  minValue<-apply(deltas,MARGIN=1,FUN=min)
  
  maxMinDiff<-maxValue-minValue
  
  precipOut$adjPrecip<<-ifelse(maxMinDiff>deltaThreshold,0,modifiedPrecip)
  
  precipOut$gaugeNoiseQF<<-ifelse(maxMinDiff>deltaThreshold,1,0)
}


##If wire deltas deviate too much from one another set the modified precip value to zero otherwise carry the modified precip value through

zeroIfDivergent(deltas=precipOut[,c("delta1","delta2","delta3")],modifiedPrecip=precipOut$intPrecip,deltaThreshold)


############################################
#USCRN code chunk - zeroIfNonpositiveDelta
#############################################

##Zero out any precip where nonpositive deltas were present from one or more of the strain gauges

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



zeroIfNonpositiveDelta(deltas=precipOut[,c("delta1","delta2","delta3")],adjustedPrecip=precipOut$adjPrecip)



###############
#Rounding data
###############

##Rounding percipitation values to the hundredth decimal place and numbers ending in 5 will be rounded up.

round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}


##Round the precipitation values
precipOut$roundedPrecip<-round2(precipOut$radjPrecip,n=2)



###############################################
# In place of USCRN code chunk - quantizeValues
###############################################

## Force any value less than 0.2 to 0
precipOut$detectablePrecip<-ifelse(precipOut$roundedPrecip<0.25,0,precipOut$roundedPrecip)

#set any depth that has an overflow flag to zero
precipOut$finalPrecip<-ifelse(precipOut$overflowQF>0,0,precipOut$detectablePrecip)



write.table(precipOut,paste("C:\\Users\\dsmith\\Desktop\\","final_priPrecip_goldenData_20160225.csv",sep="\\"),row.names=FALSE,sep=",")



#C:\\Users\\dsmith\\Desktop\\
#C:\\Users\\dsmith\\Documents\\FIU\\Git\\Met-Info\\pri_precip\\Golden Data















