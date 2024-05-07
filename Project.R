#=============================================================================
#
#           Numerical Approaches to Earth Science Problems
#
#                          Synthetic radargrams
#
#                             A.SALAMANCA
#                         GEOLOGY DEPARTMENT
#                       UNIVERSITY OF KANSAS
#                             MARCH, 2024
#
#=============================================================================

#This code generate a syntectic radargram from a fracture model with increasing aperture in R+

setwd("C:/Users/afsal/OneDrive/Desktop/KU/NumericalMethods")

#Introduce the frequency 
ff<-2e8 #200 MHz

#Introduce the time sampling rate
dtt<-2e-11  #0.2 ns/sample

#Introduce the number of samples per scan
ntt<-6000

#Introduce the velocity to make the conversion from time to depth 
velocityinmeters_s<-1.3e8

fnyq<-1.0/(2.0*dtt) #Nyquist Frequency
df<-1.0/(ntt*dtt) # Determine Fundamental Frequency

#Model of relative electrical permittivities (Model of the Earth with an heterogeneity)
ModelEarthRelativePermeabilities<-read.csv(file="ModelEarthTracesHeterogeneous.csv")

#Using a function from other file
source("Rickerwavelet.R")

impulse<-as.numeric(rickerwavelet(dtt,ntt,ff)) #Generate the time series of a Ricker wavelet

#From a Model of relative electrical permittivities, we calculate the distribution of reflection coefficients
CRS<-NA #Initialize the variable CRS (Reflection coefficients)
CRS<-rep(0, times=2*(length(ModelEarthRelativePermeabilities[,1]))) #Generates a vector of zeros of twice the length of the ricker time series
CRS<-as.data.frame(CRS) #Convert CRS to data frame

#The following loop calculates the reflection coefficient series
for (i in 1:(length(ModelEarthRelativePermeabilities[1,]))){
  for (n in 2:(length(ModelEarthRelativePermeabilities[,1]-1))){
    CRS[n,i]<-0
    CRS[n+(length(ModelEarthRelativePermeabilities[,1])),i]<-(sqrt(ModelEarthRelativePermeabilities[n+1, i])-sqrt(ModelEarthRelativePermeabilities[n, i]))/(sqrt(ModelEarthRelativePermeabilities[n+1,i])+sqrt(ModelEarthRelativePermeabilities[n,i]))
    }
  CRS[(length(ModelEarthRelativePermeabilities[,1]))*2,i]<-0
  CRS[1,i]<-0
  CRS[(length(ModelEarthRelativePermeabilities[,1])+1),i]<-0
}

#Convolution
#Prepare a data set to storage the values of the convolution serie

serie1<-NA #Initialize the matrix of the radargram
serie1<-as.data.frame(serie1)

#Add zeros to the ending of the wavelet to allow the dot product calculation in R
for (n in 1:(length(ModelEarthRelativePermeabilities[,1]))){
  impulse<-append(impulse,0)
}
impulse<-as.data.frame(impulse) #Convert the wavelet to a data frame

#Add zeros to the beginning of the reflectivity coefficient matrix to allow the dot product 
#calculation in R. The idea is having zeros tails at the beginning and ending of the series 
#to just move one position per cycle and calculate the corresponding dot product between the series, keeping the same size
#Calculate dot product between wavelet and the reflectivity coefficient series 
for (trace in 1:(length(ModelEarthRelativePermeabilities[1,]))){
  movingrow<-CRS[,trace]
  for (sample in 1:(2*(length(ModelEarthRelativePermeabilities[,1])))){
  movingrow<-append(movingrow,0) 
  movingrow<-movingrow[-1]
  serie1[sample,trace]<-sum(impulse*movingrow)
  }
}
#Generate an index from 1 to the length of the time series to plot as x value 
xindex<-0
xindex<-as.data.frame(xindex)
for (i in 1:(2*(length(ModelEarthRelativePermeabilities[,1])))){
  xindex[i,1]<-i
}
x<-as.vector(xindex[,1])#Convert as a vector the x index value
xindepth <- (x*(2e-11)*velocityinmeters_s)/2

#Plot traces, with not axes and values. Just the time series, vertically.Like a wiggle radargram 
par( mfrow= c(1,42),mar=c(0,0,0,0))

#Creating an axis to the left of the radargram in terms of time in ns
plot((as.vector(serie1[,1])),x/50, type= "n", ylab="time (ns)",ylim = rev(range(x/50)),yaxt = "n")
axis(4, pos = -1, outer = TRUE) 

for (Seriescolumn in 1:(length(ModelEarthRelativePermeabilities[1,]))){
plot((as.vector(serie1[,Seriescolumn])),x, xlab = Seriescolumn, ylim = rev(range(x)),yaxt = "n")
}

#Creating an axis to the right of the radargram in terms of meters in depth
plot((as.vector(serie1[,1])),xindepth, type= "n", ylab="y-axix label",ylim = rev(range(xindepth)),yaxt = "n")
axis(4, pos = -2) 
box() 



xindex<-0
xindex<-as.data.frame(xindex)
for (i in 1:(2*(length(ModelEarthRelativePermeabilities[,1])))){
  xindex[i,1]<-i
}

#Identify the maximum values of each trace
MaximumValues<-NA
for (Seriescolumn in 1:(length(ModelEarthRelativePermeabilities[1,]))){
  MaximumValues<-append(MaximumValues,(max(serie1[,Seriescolumn])))
}
MaximumValues<-MaximumValues[-1]
MaximumValues<-as.data.frame(MaximumValues)

#Generate an index to plot maximum values vs trace number
traceindex<-0
traceindex<-as.data.frame(traceindex)
for (trace in 1:(length(ModelEarthRelativePermeabilities[1,]))){
  traceindex[trace,1]<-trace 
}
MaximumValues$Trace<-traceindex

#Using ggplot, plot the value of the maximum amplitude vs trace
library(ggplot2)
ggplot(data = MaximumValues)
ggplot(
  data = MaximumValues,
  mapping = aes(x = Trace$traceindex, y = MaximumValues)
) +
  geom_point()


###################################-----------------################################