# Supplementary Materials: Applying Bayesian spatio-temporal models to fisheries bycatch in the Canadian Arctic R codes for the final model, inferences, and predictions

### Please refer to the INLA website: www.r-inla.org for additional information 
rm(list = ls(all = TRUE))

##################################
### MAP of BAFFIN BAY (Figure 1a)
##################################
# load datasets
#Greenland shark bycatch data
data <- read.table('Data/fackdata.txt', header=TRUE) 

# Canadian eez boundary
library(rgdal)
eez <- readOGR(dsn="Data/EEZ", layer='Atleez') 

# marine protected area bounding box 
MPAlat<-c(68.25,68.25,67.25,67.25,68.25)
MPAlong<-c(-58.551306,-60.5,-60.5,-57.8425,-58.551306)
mpa<-cbind(MPAlong,MPAlat)

library(marmap)
iso<-read.bathy("Data/etopo1_bedrock.xyz",sep=" ")
### data from NOAA: http://maps.ngdc.noaa.gov/viewers/wcs-client/

library(maps)
library(mapdata)
map('worldHires',xlim=c(-68,-50),ylim=c(66,71),fill=T, col="grey")
map.axes()
plot(iso,n=10,drawlabels=T,add=T,col='grey')
plot(eez,col="black",lty=4,add=T,lwd=2) # Canada's Eclusive Economic Zone
lines(mpa)# marine protected area
points(data$longitude,data$latitude,pch=20) # fishing locations


##########################
### Zero-inflated negative binomial 
### Bayesian spatio-temporal model 
##########################
library(INLA)
inla.setOption(inla.call='remote') 

### create mesh
mesh <- inla.mesh.2d(cbind(data$X, data$Y), max.edge=c(40, 80), cut=5)

### define SPDE
spde <- inla.spde2.matern(mesh)

### Observation matrix $A$ for each year
table(repl <- data$year-2007)#year index from 1:4
dim(A <- inla.spde.make.A(mesh, repl=repl, loc=cbind(data$X, data$Y)))

### Spatial index 
ind <- inla.spde.make.index(name='s', n.spde=spde$n.spde, n.repl=4)

### Stack data
head(data,1)
# y is Greenland shark bycatch (counts)
# covariates are specified by data[,c(4:7)]

stke <- inla.stack(data=list(gr=data$y), 
                   tag='est', A=list(A, 1), 
                   effects=list(ind, 
                       data.frame(b0=1, data[,c(4:7)])))

### INLA formula 
### inla.group reduce to unique values
form <-  gr ~ 0 + b0 + I(duration) + Ngillnets + I(TC.tspp) + 
    f(inla.group(bathymetry), model='rw1') + 
    f(s, model=spde, replicate=s.repl) 

### call INLA 
rl.nb2.gr <- inla(form, family='zeroinflatednbinomial2', 
                 data=inla.stack.data(stke), 
                control.compute=list(dic=TRUE), 
                 control.inla=list(strategy='laplace'), 
                 control.predictor=list(A=inla.stack.A(stke)),verbose=TRUE)

##########################
### RESULTS 
##########################
#Extract parameters of the random field
rf.gr <- inla.spde2.result(rl.nb2.gr, 's', spde) 

### range (km) and 95% credible intervals
c(mean=inla.emarginal(function(x) x, rf.gr$marginals.range[[1]]), 
  q=inla.hpdmarginal(0.95, rf.gr$marginals.range[[1]]))[c(1,2,3)]

#Summary table of the posterior distributions of all parameters
tables <- t(sapply(c(rl.nb2.gr$marginals.fix, 
                     rl.nb2.gr$marginals.hy[1:2], 
                     rf.gr$marginals.range[1]), function(m) 
                       c(mean=inla.emarginal(function(x) x, m), 
                         lim=inla.hpdmarginal(0.95, m)))) 
round(tables, 4)

##################################
### MAP of random field mean and sd (Figure 6)
##################################
proj<-inla.mesh.projector(mesh,xlim=range(data$X), ylim=range(data$Y))
dat = as.data.frame(proj$lattice$loc)
ndat<-data.frame(sapply(dat, rep.int, times=4))
# dim(ndat)#100 x 100 replicated 4 times for each year
ndat$year<-rep(1:4,each=10000)
table(ind$s.repl)#476 spde nodes every year
length(rl.nb2.gr$summary.random$s$mean)#1904/4

xmeans <- list()
xsd <- list()
     for (j in 1:4) {
         xmeans[[j]] <- inla.mesh.project(proj,rl.nb2.gr$summary.random$s$mean[ind$s.repl==j])
        xsd[[j]] <- inla.mesh.project(proj,rl.nb2.gr$summary.random$s$sd[ind$s.repl==j])
     }
   
ndat <- data.frame(ndat, unlist(xmeans),unlist(xsd)) 
head(ndat)
ndat$year <- factor(ndat$year+2007)
colnames(ndat)[4:5]<-c("mean","sd")

### to plot bathy
r1 <- as.raster(iso)
proj <- "+proj=longlat +datum=WGS84"
library(raster)
r2 <- projectRaster(r1,crs=proj)
Arcticproj<-"+proj=utm +zone=19 +ellps=GRS80 +units=km +no_defs"
r.km <- projectRaster(r1,crs=Arcticproj)

library(fields)
library(lattice)
### for grey colour scale: col.regions=gray.colors(16,start=1,end=0)
my.at <- seq(-3500,0,500)
levelplot(mean ~  V1*V2 | year, data=ndat,xlab='', ylab='',
       col.regions=tim.colors(100),scale=list(draw=FALSE),par.strip.text=list(cex=2),strip = strip.custom(bg="white")) +  latticeExtra::layer(sp::sp.polygons(land.km,fill="lightgrey"))+rasterVis::contourplot(r.km, labels=list(cex=1), at=my.at,col="black")

levelplot(sd ~  V1*V2 | year, data=ndat,xlab='', ylab='',
          col.regions=tim.colors(100),scale=list(draw=FALSE),par.strip.text=list(cex=2),strip = strip.custom(bg="white")) +  latticeExtra::layer(sp::sp.polygons(land.km,fill="lightgrey"))+rasterVis::contourplot(r.km, labels=list(cex=1), at=my.at,col="black")

##########################
### PREDICTIONS
##########################
#load grid data
dim(bgrid <- read.table('Data/gridpred.txt', header=TRUE))

#*Grid locations* We used our final model to predict expected bycatch in areas 
#neigbouring observed fishing hauls. However, since the fishery expended 
#spatially, the spatial fishing pattern changed over time with some years with 
#no fishing in some areas of our domain (e.g., northeast cluster at 71 degrees 
#was not fished in 2008 and 2009). For these area-time, we did not predict
#bycatch. However, note that it would be possible to do so (assuming some
#average fishing conditions (duration, Ngillnets, TC.tspp), but the spatial
#random effect would be very close to zero and with large variance. We selected
#nearest neighbour grid points less than 30km from our fishing hauls
#using nndistF function from [splancs].

# par(mfrow=c(1,1))
# with(data, plot(X, Y, asp=1,pch='.'))
# points(bgrid$x, bgrid$y, col=gray(.7), pch=4, cex=0.01)

library(splancs) 
iig <- sapply(2008:2011, function(y) {
    ii.y <- which(data$year==y)
    d <- nndistF(cbind(data$X[ii.y], data$Y[ii.y]), cbind(bgrid$x, bgrid$y))
    which(d<30)
})
str(iig)

# *Prediction of Greenland halibut (TC.tspp)* Since Greenland halibut catch (kg)
# are driven by unknown variable, we also interpolated TC.tspp over the grid for
# each year. The interpolation was done using a weigthed mean, where the weights
# are proportional to exp(-dist/a), where 'a' is equal to 0.5km, that is inversely
# proportional to the total grid distance (1x1km).

library(geoR)
itc <- lapply(1:4, function(y) {
    ii.y <- which(data$year==((2008:2011)[y]))
    d <- loccoords(cbind(bgrid$x, bgrid$y)[iig[[y]],], 
                   cbind(data$X, data$Y)[ii.y,]) 
    w <- exp(-d/0.5)
    drop((w/rowSums(w))%*%(data$TC.tspp[ii.y]))
})
sapply(itc, summary)
summary(data$TC)

### Building prediction stacks
### for each year with
###  duration=15.94, Ngillnets=43.44
###  bathymetry over grid and TC.tspp interpolated
y.stk.p <- lapply(1:4, function(y) {
    cat('year', (2008:2011)[y], '\n')
    Ap <- inla.spde.make.A(mesh, cbind(
        bgrid$x[iig[[y]]], bgrid$y[iig[[y]]]), repl=y, n.repl=4)
    inla.stack(data=list(gr=NA), 
               tag=paste('prd', y, sep=''), A=list(Ap, 1), 
               effects=list(ind, 
                   data.frame(b0=1, bathymetry=bgrid$bathymetry[iig[[y]]],
                              duration=15.98, Ngillnets=43.44,
                              TC.tspp=itc[[y]])))
})

### we can join all stack data into one and run it all at once, 
### but this is expensive. So we use each prediction stack per year.
yprd <- lapply(y.stk.p, function(s) {
    stk <- inla.stack(stke, s) 
    res <- inla(form, family='zeroinflatednbinomial2', 
                data=inla.stack.data(stk), quantiles=NULL, 
                control.predictor=list(A=inla.stack.A(stk), compute=TRUE, quantiles=NULL), 
                control.compute=list(return.marginals=FALSE), 
                inla.call='remote')
    print(res$cpu)
    ii <- inla.stack.index(stk, names(s$effects$index))$data 
    res$summary.linear.pred[ii,]
})

### join all into data.frame
jprd <- Reduce('rbind', lapply(1:4, function(y) data.frame(
    x=bgrid$x[iig[[y]]], y=bgrid$y[iig[[y]]], lpred=yprd[[y]]$mean, year=(2008:2011)[y])))

# write.table(jprd, 'Data/jprd.txt', row.names=FALSE) #save predictions

##################################
### MAPs of Prediction (Figure 7)
##################################
#read the prediction outputs
jprd <- read.table('Data/jprd.txt',header=T)
jprd.y <- split(jprd, jprd$year)
str(jprd.y)

proj <- "+proj=utm +zone=19 +ellps=GRS80 +units=km +no_defs" 
res0 <- lapply(jprd.y, function(x) SpatialPoints(cbind(x$x, x$y), CRS(proj)))

newproj <- "+proj=longlat +zone=19 +ellps=GRS80" # lat/long
res <- lapply(1:length(res0), function(i)
              SpatialPointsDataFrame(spTransform(res0[[i]], CRS(newproj)),
                              data.frame(pred=exp(jprd.y[[i]]$lpred))))

q <- c(0, 1, 3, 5, 10, 30)
leg <- paste(paste(c('<', q[2:4], '>'), colapse='', sep=''),
             c('', rep('-', 3), ''),
             paste(c(q[2:5], q[5]), colapse='', sep=''), sep='')

par(mfrow=c(2,2), mar=c(2,2,0.5, 0.5), mgp=c(2,1,0))
for (i in 1:4) {
    if (i != 3){
    map('worldHires',xlim=c(-68,-50),ylim=c(66,71),col="lightgrey")
    box()
    plot(iso,n=10,drawlabels=T,add=T,col='grey')
    plot(eez,col="black",lty=4,add=T,lwd=2) 
    lines(mpa)
    ctest<-tim.colors(6)
    plot(res[[i]], add=T, col=ctest[findInterval(res[[i]]$pred, q)], cex=0.1, pch=19)
    text(-65,66.5,labels=eval(i+2007),font=2,cex=1.2)
  }
  else {
    map('worldHires',xlim=c(-68,-50),ylim=c(66,71),col="lightgrey")
    map.axes()
    plot(iso,n=10,drawlabels=T,add=T,col='grey')
    plot(eez,col="black",lty=4,add=T,lwd=2) 
    lines(mpa)
    plot(res[[i]], add=T, col=ctest[findInterval(res[[i]]$pred, q)], cex=0.1, pch=19)
    text(-65,66.5,labels=eval(i+2007),font=2,cex=1.2)
    legend(-57, 71, leg, fill=ctest,bg="white")
  }
}

##################################
### MAP of BAFFIN BAY + MESH (Figure 1b)
##################################
# Project to km as mesh
land<-map("worldHires", xlim=c(-68,-50),ylim=c(66,71), fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(land$names, ":"), function(x) x[1])
library(maptools)
land.sp <- map2SpatialPolygons(land, IDs=IDs, proj4string=CRS("+proj=longlat +datum=WGS84"))

land.km <- spTransform(land.sp, CRS(Arcticproj))

plot(mesh,asp=1,main=NULL,sub=NULL)
plot(land.km,add=T)
box(lwd=2)

