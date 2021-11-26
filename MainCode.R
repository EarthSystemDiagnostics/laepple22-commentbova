
## Code to reproduce the results/figures in Laepple et al., Nature 2022; Concerns of assuming linearity in the reconstruction of thermal maxima;
## Note that some formatting details as the axis labeling etc; are not fully
##reproduced as a tradeoff of keeping a shorter code
## The code is not optimized for efficiency / speed; Figure 1 should take <10min; Figure 2 < 1h.
## For questions, contact tlaepple@awi.de 



library(ncdf4)
library(RColorBrewer)
FILEPATH_SST<-"/Users/tlaepple/data/bova/data/SST_Regridded/"
FILEPATH_LIB<-"/Users/tlaepple/data/bova/src/lib/"

# install.packages("remotes")  (please use an up to date version of remotes!) 
remotes::install_github("EarthSystemDiagnostics/pfields") ##Functions to work with 3D Fields
library(pfields)
source(paste(FILEPATH_LIB,"selspace.pField.R",sep="")) ## Function to select a region in a 3D field


source(paste(FILEPATH_LIB,"mon_insolation.R",sep="")) ##Functions for insolation calculation
source(paste(FILEPATH_LIB,"satmethodfunctions.R",sep="")) ##Functions implementing the SAT method
source(paste(FILEPATH_LIB,"read_ccsmnc.R",sep="")) ##reading function for netcdf
source(paste(FILEPATH_LIB,"metrics.R",sep="")) ##Helper functions for metrics

source(paste(FILEPATH_LIB,"plotlib.R",sep="")) ##Helper functions to plot a field as map 


orbit91<-read.table(paste(FILEPATH_LIB,"ins_data.txt",sep=""))  ##Read orbital parameters and precalculate them
orbital_global<-orbital_parameters((0:50000)/10)



### Read all the climate model data provided by Feng He (stored at https://osf.io/x695d/)
eem.ann<-read_ccsmnc(paste(FILEPATH_SST,"b30t31_trs_ogg.pop_ANN_means.SST.128ka-115ka.TEMP.remap.ESMF.T31.nc",sep=""))

eem.mon<-list()
for (iMon in 1:9) {
    FILENAME<-paste(FILEPATH_SST,"b30t31_trs_ogg.pop.h.SST.0",iMon,".128ka-115ka.TEMP.remap.ESMF.T31.nc",sep="")
    eem.mon[[iMon]]<-read_ccsmnc(FILENAME)
}

for (iMon in 10:12) {
    FILENAME<-paste(FILEPATH_SST,"b30t31_trs_ogg.pop.h.SST.",iMon,".128ka-115ka.TEMP.remap.ESMF.T31.nc",sep="")
    eem.mon[[iMon]]<-read_ccsmnc(FILENAME)
}


hol.ann<-read_ccsmnc(paste(FILEPATH_SST,"b30t31_trs_ogg.pop_ANN_means.SST.12ka-0ka.TEMP.remap.ESMF.T31.nc",sep=""))


hol.mon<-list()
for (iMon in 1:9) {
    FILENAME<-paste(FILEPATH_SST,"b30t31_trs_ogg.pop.h.SST.0",iMon,".12ka-0ka.TEMP.remap.ESMF.T31.nc",sep="")
    hol.mon[[iMon]]<-read_ccsmnc(FILENAME)[2:120,] ##Cut to the same time interval as hol.ann
}

for (iMon in 10:12) {
    FILENAME<-paste(FILEPATH_SST,"b30t31_trs_ogg.pop.h.SST.",iMon,".12ka-0ka.TEMP.remap.ESMF.T31.nc",sep="")
    hol.mon[[iMon]]<-read_ccsmnc(FILENAME)[2:120,]  ##Cut to the same time interval as hol.ann
} 






### First experiment; Test if the SAT method can identify seasonal biases by applying it to the annual mean Eemian time-series
annualfit<-fitField(eem.ann)

quartz(width=4.7*1.3,height=3*1.3)
plotmap.square(annualfit[1,],shift=TRUE,main="Detected season; 'Truth' AOGCM MAT ")

sum(annualfit[1,]==13,na.rm=TRUE)/sum(!is.na(annualfit[1,]))  ##0.21 = 79% of the gridboxes are wrongly identified as seasonally biased 


annualfit.tropics <- selspace(annualfit,lat1=-40,lat2=40,lon=0,lon2=360)
sum(annualfit.tropics[1,]==13,na.rm=TRUE)/sum(!is.na(annualfit.tropics[1,]))  ##0.21... same result for the tropics

(sum(annualfit.tropics[1,]<3,na.rm=TRUE)+sum(annualfit.tropics[1,]>10,na.rm=TRUE))/sum(!is.na(annualfit.tropics[1,])) #63% Jan/Feb & Nov/Dec/Annual = 37% in March-October


## Second experiment: Test if the SAT method is skillful in correcting seasonally biased records
## Run it on the EEM model data for any of the 12 recording monts
fitResult<-list()
for (iMon in 1:12) fitResult[[iMon]]<-fitField(eem.mon[[iMon]])

corrected<-list()
for (iMon in 1:12) corrected[[iMon]]<-correctField(fitResult[[iMon]],eem.mon[[iMon]],hol.mon[[iMon]])


#### Test how often the correct proxy season (correct (defined here as two to zero months prior to the prescribed month to allow for a delay between insolation and temperature due to heat-capacity) is detected
## Global
fraction.correct<-vector()
for (iMon in 1:12) fraction.correct[iMon]<-distance_month(c(fitResult[[iMon]][1,]),iMon)/sum(!is.na(fitResult[[iMon]][1,]))

mean(fraction.correct) #0.427




##Correlate every on the 12 fields with the annual field, once before correction = cor.raw and once after correction cor.corrected
cor.corrected <- lapply(corrected,function(x) return(corField(x$hol,hol.ann)))
cor.raw <- lapply(hol.mon,function(x) return(corField(x,hol.ann)))


### Count how often the correction worked; here defined as an improvement in correlation
count <- pField(0,1,getlat(hol.ann),getlon(hol.ann),"#month where correction works")
for (iMon in 1:12)
    count[1,]<-count[1,]+(cor.corrected[[iMon]]>cor.raw[[iMon]])


## Document how often all months, now months and less than 6 monts were skillful                                     
NBoxTotal <- sum(!is.na(count))
sum(count[1,]==0,na.rm=TRUE)/NBoxTotal  #0.258
sum(count[1,]==12,na.rm=TRUE)/NBoxTotal #0.0042
sum(count[1,]<6,na.rm=TRUE)/NBoxTotal #0.74


## Repeat for 40S-40N
count.tropics <- selspace(count,lat1=-40,lat2=40,lon1=0,lon2=360)                              
NBoxTotal.tropics <- sum(!is.na(count.tropics))
sum(count.tropics[1,]==0,na.rm=TRUE)/NBoxTotal.tropics  #0.136
sum(count.tropics[1,]==12,na.rm=TRUE)/NBoxTotal.tropics #0.080
sum(count.tropics[1,]<6,na.rm=TRUE)/NBoxTotal.tropics #0.67


## Plot Figure 1
pal=colorRampPalette(rev(brewer.pal(9,"Spectral")[2:9]))

quartz(width=4.7*1.3,height=3*1.3)
pdf(width=4.7*1.3,height=3*1.3,file="./plots/Figure1_raw.pdf")
plotmap.square(count,zlim=c(0,12),shift=TRUE,stype=2,palette=pal,main="",levels=-1:12+0.5)
dev.off()

#### 


### Third experiment:


## Create a set of artificial trends
N<-dim(eem.ann)[1]
tstrend0<-seq(from=0,to=0,length.out=N) #No trend
tstrend1<-seq(from=3,to=0,length.out=N) #linear trend
tstrend2<-(cos(1:N*2*pi/N/4-0.2))*3 # piece of a cosine function
tstrend3<-(sin(1:N*2*pi/N/2))*2 # another piece of a sine/cosine function


##Create 3D fields with the different trends defined above
trend1 <- eem.ann  ## Copy eem.ann and overwrite it 
trend1[]<-rep(tstrend1,by=dim(eem.ann)[2])

trend2 <- eem.ann  ## Copy eem.ann and overwrite it 
trend2[]<-rep(tstrend2,by=dim(eem.ann)[2])

trend3 <- eem.ann ## Copy eem.ann and overwrite it 
trend3[]<-rep(tstrend3,by=dim(eem.ann)[2])


##For simplicity, the results are all stored in the same list in the index 1:12,13:24,25:36
## Fit the insolation model
fitResult.trend<-list()
for (iMon in 1:12) fitResult.trend[[iMon]]<-fitField(eem.mon[[iMon]]+trend1)
for (iMon in 13:24) fitResult.trend[[iMon]]<-fitField(eem.mon[[iMon-12]]+trend2)
for (iMon in 25:36) fitResult.trend[[iMon]]<-fitField(eem.mon[[iMon-24]]+trend3)


##and use the fit to correct the Eemian (and the Holocene, not further used here) using the SAT method; the eemian input fields are the CCSM simulated EEM + the trends; the Holocene input fields are the uncorrected Holocene simulations; 
corrected.trend<-list()
for (iMon in 1:12) corrected.trend[[iMon]]<-correctField(fitResult.trend[[iMon]],eem.mon[[iMon]]+trend1,hol.mon[[iMon]])
for (iMon in 13:24) corrected.trend[[iMon]]<-correctField(fitResult.trend[[iMon]],eem.mon[[iMon-12]]+trend2,hol.mon[[iMon-12]])
for (iMon in 25:36) corrected.trend[[iMon]]<-correctField(fitResult.trend[[iMon]],eem.mon[[iMon-24]]+trend3,hol.mon[[iMon-24]])


### Figure 2
colors=colorRampPalette(rev(brewer.pal(9,"Spectral")[2:9]))(12)

lwd=2
quartz(width=6,height=2.5)
par(mfrow=c(1,4))
ylim=c(-1.6,1.6)
par(mai=c(0.6,0,0.6,0))

plot(scale(areasubmean(eem.ann),scale=FALSE),ylim=ylim,lwd=lwd,main="",axes=FALSE,xlab="",ylab="")
box()
axis(1,at=c(-126,-121,-116))

for (iMon in 1:12) lines(scale(areasubmean(corrected[[iMon]]$eem),scale=FALSE),col=colors[iMon],lwd=1)

lines(c(time(eem.ann)),scale(tstrend0,scale=FALSE),col="grey",lwd=2,lty=2)


plot(scale(areasubmean(eem.ann+trend1),scale=FALSE),ylim=ylim,lwd=lwd,main="",axes=FALSE,xlab="",ylab="")
box()
axis(1,at=c(-126,-121,-116))
for (iMon in 1:12) lines(scale(areasubmean(corrected.trend[[iMon]]$eem),scale=FALSE),col=colors[iMon],lwd=1)
lines(c(time(eem.ann)),scale(tstrend1,scale=FALSE),col="grey",lwd=2,lty=2)

plot(scale(areasubmean(eem.ann+trend2),scale=FALSE),ylim=ylim,lwd=lwd,main="",axes=FALSE,xlab="",ylab="")
box()
axis(1,at=c(-126,-121,-116))
for (iMon in 13:24) lines(scale(areasubmean(corrected.trend[[iMon]]$eem),scale=FALSE),col=colors[iMon-12],lwd=1)

lines(c(time(eem.ann)),scale(tstrend2,scale=FALSE),col="grey",lwd=2,lty=2)


plot(scale(areasubmean(eem.ann+trend3),scale=FALSE),ylim=ylim,lwd=lwd,main="",axes=FALSE,xlab="",ylab="")
box()
axis(1,at=c(-126,-121,-116))
for (iMon in 25:36) lines(scale(areasubmean(corrected.trend[[iMon]]$eem),scale=FALSE),col=colors[iMon-24],lwd=1)

lines(c(time(eem.ann)),scale(tstrend3,scale=FALSE),col="grey",lwd=2,lty=2)

###




     

