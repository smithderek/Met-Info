
#Thresholds
#The largest step in precip allowed between time stamps
#25

#The largest difference allowed among the wire deltas
#0.5

## Low depth range
#-10


######################################################################################
## Simulate precipitation data for three hour periods to be used ot test the algorithm
######################################################################################


##Number of precip values to generate
size<-5000

#baseline precip level
baseline<-42

#Rain intensity
#low
lam<-0.001

#Med
lam<-0.005

#High
lam<-0.01

#generate precip frequency
set.seed(1234)
precip_freq<-rpois(n=size, lambda=lam)

#generate precip amounts
set.seed(1234)
precip_amt<-rnorm(size,mean=0.5,sd=0.2)

#expected preip
expPrecip<-precip_amt*(1.5*precip_freq)
plot(expPrecip)
#actual measured precip based on 0.2 mm threshold
measPrecip<-ifelse(expPrecip> 0.2,expPrecip,0)

#cumulative sum of precips
cumPrecip<-cumsum(expPrecip)

#add on initial starting point for gauge
precipTotals<-cumPrecip+baseline

#plot precip totals
plot(precipTotals)
summary(precipTotals)

#####################
####Noise level
####################

########################################
##None
noiseLow1<-0
noiseLow2<-0
noiseLow3<-0


#low
#gauge 1
set.seed(1234)
noiseLow1<-rnorm(size,mean=0,sd=0.01)
#gauge 2
set.seed(2345)
noiseLow2<-rnorm(size,mean=0,sd=0.01)
#gauge 3
set.seed(3456)
noiseLow3<-rnorm(size,mean=0,sd=0.01)
########################################

#Med
set.seed(1234)
noiseMed<-rnorm(size,mean=0,sd=0.1)
#High
set.seed(1234)
noiseHi<-rnorm(size,mean=0,sd=0.5)

noiseM<-rnorm(size,mean=0,sd=0.1)
#add noise to measurement
lowNoisePrecip<-precipTotals+noiseLow
plot(lowNoisePrecip)

medNoisePrecip<-precipTotals+noiseMed
plot(medNoisePrecip)


highNoisePrecip<-precipTotals+noiseHi
plot(highNoisePrecip)

##############################################
### Add noise to wires
#############################################
wire1<-noiseLow1+precipTotals

wire2<-noiseLow2+precipTotals

wire3<-noiseLow3+precipTotals



###########################################
## Convert to frequency values
#########################################
#
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


######################
##  Determine frequency from depth
######################
a1<-c2_1
b1<-(-2*c2_1*f0_1)+c1_1


a2<-c2_2
b2<-(-2*c2_2*f0_2)+c1_2

a3<-c2_3
b3<-(-2*c2_3*f0_3)+c1_3

c1<-(-f0_1*c1_1)+(c2_1*f0_1^2)-wire1
wire1Freq<-((-1*b1)+sqrt(b1^2-4*a1*c1))/(2*a1)

c2<-(-f0_2*c1_2)+(c2_2*f0_2^2)-wire2
wire2Freq<-((-1*b2)+sqrt(b2^2-4*a2*c2))/(2*a2)


c3<-(-f0_3*c1_3)+(c2_3*f0_3^2)-wire3
wire3Freq<-((-1*b3)+sqrt(b3^2-4*a3*c3))/(2*a3)





wireFreq<-data.frame(cbind(wire1Freq,wire2Freq,wire3Freq))



##############################################
## REPLICATe values for a five minute period
##############################################

repFreq<-wireFreq[rep(row.names(wireFreq),times=30),]

##create a numeric vector of row numbers to sort by
repFreq$rowNum<-as.numeric(row.names(repFreq))

##Sort the data set by the row numbers
finalFreq<-repFreq[order(repFreq$rowNum),]


################################
## ADD in Jitter TBD
################################

#jfin<-jitter(wireFreq[,1])

#####################################

#################################
## Add timestamps frequencies
#################################
# make date and time



start <- as.POSIXct("2012-01-15",tz="UTC")
interval <- 10

##Time formatted for CI
time<-format(seq(from=start, by=interval, length.out=nrow(finalFreq)),"%d-%b-%y %I.%M.%OS3 %p")

##Normal R time format
#time<-seq(from=start, by=interval, length.out=nrow(finalFreq))


precip<-data.frame(time=time,Freq1=finalFreq$wire1Freq,Freq2=finalFreq$wire2Freq,Freq3=finalFreq$wire3Freq)

#colnames(precip)<-c("MEAS_RDOT_STRT_DATE","N_VAL","VAL_ID")

write.table(precip[,c(1,4)],"C:\\Users\\dsmith\\Documents\\FIU\\Precipitation\\Primary Precip Gauge\\USCRN data\\Golden Data\\priPrecip_clean_goldenData_maxID_13824_gauge3.csv",row.names=FALSE,sep=",")


##create one minute means over the time series
precip.min.means <- aggregate(precip[,c(2:4)], list(time=cut(precip$time,breaks="min")),mean,na.rm = TRUE)

## convert time stamps to posixct
precip.min.means$time<-as.POSIXct(precip.min.means$time)

##Subset the one minute means at five minute increments

#create a vector for the minutes
precip.min.means$min<-format(precip.min.means$time,"%M")


#subset the data set at five minute increments
x<-precip.min.means[precip.min.means$min=="00" |precip.min.means$min=="05"|precip.min.means$min=="10"|precip.min.means$min=="15"|precip.min.means$min=="20"|precip.min.means$min=="25"|precip.min.means$min=="30"|precip.min.means$min=="35"|precip.min.means$min=="40"|precip.min.means$min=="45"|precip.min.means$min=="50"|precip.min.means$min=="55",]


##################################################
## Start adding in algorithm to convert to depths
#################################################
head(x)

wire1D<-c1_1*(x$Freq1-f0_1)+(c2_1*((x$Freq1-f0_1)^2))
wire2D<-c1_2*(x$Freq2-f0_2)+(c2_2*((x$Freq2-f0_2)^2))
wire3D<-c1_3*(x$Freq3-f0_3)+(c2_3*((x$Freq3-f0_3)^2))

out<-data.frame(wire1D,wire2D,wire3D)


precip<-read.table("C:\\Users\\dsmith\\Documents\\FIU\\Precipitation\\Primary Precip Gauge\\USCRN data\\Golden Data\\priPrecip_clean_goldenData_maxID_13824_all.csv",header=TRUE,sep=",")

pre<-precip[,2]

######################################
## Add gaps to data set
#######################################


gap<-function(data,nGaps=1,minGap=1,maxGap=1){
  
  ##vector the length of the data set
  len<-1:length(data)
  
  #length of data set
  end<-length(data)
  
  ##gap indicator
  gaps<-rep(0,length(data))
  
  #data with gaps
  data.gap<-data
  
  for(i in 1:nGaps){
    dur<-sample(minGap:maxGap,size=1) #duration of gap
    gapStart<-sample(len, size=1) #Start of gap in data
    gapEnd<-gapStart+dur #End point of gap in data
    
    if(gapEnd <= end){
      data.gap[gapStart:gapEnd]<-NA #input NA into data set for the gap
      gaps[gapStart:gapEnd]<-1 #update indicator variable for the gap period
      
    }
    
    if(gapEnd > end){
      gapEnd<-end
      data.gap[gapStart:gapEnd]<-NA #input NA into data set for the gap
      gaps[gapStart:gapEnd]<-1 #update indicator variable for the gap period
      
    }
    
    
  }
  
 out<<-data.frame(data=data,data.gap=data.gap,gaps=gaps)
  
  
}


gap(data=pre,nGaps=5,minGap=1,maxGap=20)

out$time<-precip$time

h2 <- ggplot(data=out,aes(y=data,color=gaps))+ 
  geom_point(size=2.5)+
  #scale_color_manual(values = c('blue', 'red'))+ 
  scale_colour_gradient(limits=c(0,65),low="blue", high="red",name="Temp Deg C")+
  geom_errorbar(h2limits, width=0.2,color="black")+
  ylab("heated-nonheated (ppt)")+
  xlab("Time")+
  ggtitle("Deuterium- 5 minute averages durring heating events")
h2


h2 <- ggplot(data=out,aes(y=data,x=time))+
  geom_line(size=2.5)
h2







