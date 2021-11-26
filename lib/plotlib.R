#some plotting routines for fields

plotsquare<-function (x = seq(0, 1, len = nrow(z)), y = seq(0, 1, len = ncol(z)), 
    z, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), 
    zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels), 
    nlevels = 20, color.palette = cm.colors, col = color.palette(length(levels) - 
        1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
    xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, 
    ...) 
{
    if (missing(z)) {
        if (!missing(x)) {
            if (is.list(x)) {
                z <- x$z
                y <- x$y
                x <- x$x
            }
            else {
                z <- x
                x <- seq(0, 1, len = nrow(z))
            }
        }
        else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
        y <- x$y
        x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
        stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    layout(matrix(c(2, 1), nc = 2), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4] <- mar[2]
    mar[2] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
        yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1], col = col)
    if (missing(key.axes)) {
        if (axes) 
            axis(4)
    }
    else key.axes
    box()
    if (!missing(key.title)) 
        key.title
    mar <- mar.orig
    mar[4] <- 1
    par(mar = mar)
    #plot.new()
    #plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    image(x,y,z,col=col)
    if (missing(plot.axes)) {
        if (axes) {
            title(main = "", xlab = "", ylab = "")
            axis(1)
            axis(2)
        }
    }
    else plot.axes
    if (frame.plot) 
        box()
    if (missing(plot.title)) 
        title(...)
    else plot.title
    invisible()
}


myfun1<-function(sTitle,sSub)
{
       title(main=sTitle,sub=sSub)
       addland(col="black")
       grid()
}

myfun2<-function(sTitle,sSub,lat,lon,plotdata)
{ title(main=sTitle,sub=sSub)
contour(lon,lat,plotdata,col="white",zlim=c(min(plotdata),0),lty=1,lwd=2,nlevels=5,add=TRUE);
contour(lon,lat,plotdata,col="grey20",zlim=c(0,max(plotdata)),lty=1,lwd=2,nlevels=5,add=TRUE);
addland(col="black")
grid()
}

#checks the data
#wrap the data to get a continous lat/lon field
#reverse latitudes if needed
plot.preparation <- function(plotdata)
{
        temp<-attributes(plotdata)
  	if (prod(dim(plotdata)) != length(temp$lon)*length(temp$lat)) stop("N(data) != N(lat)*N(lon)")

        plotdata<-matrix(plotdata,length(temp$lon),length(temp$lat)) #make a 2D array
        
        #arrange to get a continous field
        d<-diff(temp$lon)
        if (max(d) > (min(d)+0.01))
          { nlon<-length(temp$lon)
            edgelon<-which(d==max(d))
            plotdata<-rbind(plotdata[(edgelon+1):nlon,],plotdata[1:edgelon,])
            temp$lon<-c(temp$lon[(edgelon+1):nlon],temp$lon[1:edgelon]+360)
          }
        
        
        if (temp$lat[2] < temp$lat[1])    #if the latitudes are from + to -, reverse them
         {
          temp$lat<-rev(temp$lat)
          plotdata<-plotdata[,rev(seq(len=ncol(plotdata)))]
        }
        return(list(data=plotdata,lat=temp$lat,lon=temp$lon))
}
rbow <- function (n, s = 1, v = 1, start = 0, end = 0.7, gamma = 1) 
{
    if ((n <- as.integer(n[1])) > 0) {
        if (start == end || any(c(start, end) < 0) || any(c(start, 
            end) > 1)) 
            stop("'start' and 'end' must be distinct and in [0, 1].")
        result<-hsv(h = seq(start, ifelse(start > end, 1, 0) + end, length = n)%%1, 
            s, v, gamma)
        rev(result)    
    }
    else character(0)
}

plotmap.square <- function(plotdata,main=NULL,zlim=range(plotdata,finite=TRUE),levels=pretty(zlim,nlevels),nlevels=20,palette=rbow,FUN=NULL, ...)
{
	temp<-attributes(plotdata)
      sSub<-NULL
	if (time(plotdata) != 9999) sSub<-paste("time:",format(time(plotdata)))
	if (is.null(main)) main<-temp$name
        tmp<-plot.preparation(plotdata)
	plotsquare(tmp$lon,tmp$lat,tmp$data,zlim=zlim,nlevels=nlevels,levels=levels,color=palette,plot.title={
        title(main=main,sub=sSub);
        if (!is.null(FUN)) FUN(tmp$lon,tmp$lat,tmp$data)
        addland(col="black");
        grid()})
      }


