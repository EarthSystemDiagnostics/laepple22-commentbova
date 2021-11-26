
##' @title Correlate two time-varying fields by correlating each gridbox in field 1 with the corresponding box in field 2
##' @param field1 pField object: single 3D field
##' @param field2 pField object: single 3D field
##' @return field pField object with the correlation coefficients
##' @author Thomas Laepple 
corField<-function(field1,field2)
{
    index <- (!is.na(field1)[1,])
    result<-pField(NA,1,getlat(field1),getlon(field1),"Correlation")
    for (i in which(index)) {
        if (i%%100) print(paste(i,",",length(index)))
        result[,i]<-cor(field1[,i],field2[,i])
    }
    return(result)
    
}



## Distance measure
##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param time1 
##' @param time2 
##' @return 
##' @author Thomas Laepple 
distance_month <- function(time1,time2)
{
    distance <- (time2-time1)
    index <- which(distance<0)
    distance[index]<-distance[index]+12
    return(sum(distance<=2,na.rm=TRUE))
    
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
##' @param field 
##' @param lat1 
##' @param lat2 
##' @return 
##' @author Thomas Laepple 
areasubmean<-function(field,lat1=-40,lat2=40)   areamean(selspace(field,lat1=lat1,lat2=lat2,lon1=0,lon2=360))


