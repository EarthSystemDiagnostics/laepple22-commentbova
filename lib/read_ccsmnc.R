
##' @Reading function for the CCSM files
##' @param FILENAME FILENAME of the CCSM netcdf file
##' @param varname name of the variable in the netcdf file
##' @param name of the variable in the pField object
##' @return pField object 
##' @author Thomas Laepple 
read_ccsmnc<-function(FILENAME="",varname="TEMP",name="")
{
  temp.nc = nc_open(FILENAME)
  #Read out the data
  temp.time <- ncvar_get(temp.nc,"time")
  temp.date<-temp.time
  temp.data <-ncvar_get(temp.nc,varname)
  temp.lat <-ncvar_get(temp.nc,"lat")
  temp.lon <-ncvar_get(temp.nc,"lon")
  #Sort the latitudes
  tmp<-sort(temp.lat,index.return=TRUE)
  temp.lat<-temp.lat[tmp$ix]
  temp.data<-temp.data[,tmp$ix,]
  return(pField(temp.data,temp.date,lat=temp.lat,lon=temp.lon,name=name,history=FILENAME))
}
