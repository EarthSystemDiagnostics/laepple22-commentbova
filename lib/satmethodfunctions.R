
## Functions implementing the Bova et al., SAT method for time-series and fields

##' @title Fit the time-series to any monthly and annual insolation
##' @param tsProxy (pseudo)proxy time series as pTs object (with time in kyr BP, and latitude in degN)
##' @return vector with (best fitting month, slope in K/(Wm^-1), correlation coef)
##' ## If the annual mean is fitting best, 13 is returned as best fitting month
##' @author Thomas Laepple 
fit_insolation<-function(tsProxy)
{
    insol<-mon_insolation(-1*time(tsProxy),getlat(tsProxy))
    ## First 'detect seasonal bias' by maximising the correlation 
    corvals <- cor(tsProxy,insol)
    month<-which.max(corvals)
    if (month == 13) return(c(13,0,corvals[month]))  ##Return 0 slope if the best fit is annual

    ## If there is a seasonal bias (month != 13) than estimate and store the slope to the insolation anomaly
    delta.insol <- (insol[,month]-insol[,13])
    model<-lm(c(tsProxy)~delta.insol)
    slope<-model$coeff[2]  ##called As in Bova et al.,
    return(c(month,slope,corvals[month]))
}



##' @title 'Correct' timeseries using the SAT method
##' @param pTs object of the time-series to be corrected 
##' @param fit vector with three elements  (best fitting month, slope in K/(Wm^-1), correlation coef)
##' @return pTs object: corrected timeseries
##' @author Thomas Laepple 
correct_record <- function(tsProxy,fit)
{
    insol<-mon_insolation(-1*time(tsProxy),getlat(tsProxy))
    delta.insol <- insol[,fit[1]]-insol[,13] ##Called deltaI in Bova
    orbitalBias<-fit[2]*delta.insol  ##Called deltaSST in Bova
    return(tsProxy-orbitalBias)
}





##' @title Fit every time-series in a field to the insolation 
##' @param climfield pField object
##' @return pField object with three fields: best fitting month, slope, correlation to the insolation at the best fitting month
##' @author Thomas Laepple 
fitField <- function(climfield)
{
    fitResult <- pField(NA,1:3,lat=getlat(climfield),lon=getlon(climfield),name=c("BestMonth,Slope,R"))
    index <- (!is.na(climfield)[1,]) ##Only apply it to non-missing = ocean data
    for (i in which(index)) {
        if (i%%100) print(paste(i,",",length(index))) 
        tsClimate <- climfield[,i]
        tsProxy<-tsClimate
        fitResult[,i]<-fit_insolation(tsProxy)
    }
    return(fitResult)
}


    


##' @title Apply the SAT method on the Eemian and the Holocene model fields 
##' @param fitField pField object with output of the fitting
##' @param clim.eem pField object with the Eemian climate 
##' @param clim.hol pField object with the Holocene climate 
##' @return list with hol = corrected Holocene climate, eem = corrected Eemian climate
##' @author Thomas Laepple 
correctField <- function(fitField,clim.eem,clim.hol)
{
    ##Build empty fields with the same dimensions as the uncorrected ones.
    eem.corrected <- clim.eem
    eem.corrected[]<-NA

    hol.corrected <- clim.hol
    hol.corrected[]<-NA
    
    index <- (!is.na(clim.eem)[1,])
    for (i in which(index)) {
        
        if (i%%100) print(paste(i,",",length(index)))
        
        eem.corrected[,i]<-correct_record(clim.eem[,i],c(fitField[,i]))
        hol.corrected[,i]<-correct_record(clim.hol[,i],c(fitField[,i]))

    }
    return(list(eem=eem.corrected,hol=hol.corrected))
}

