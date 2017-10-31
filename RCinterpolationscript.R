###Initial script for RC project
###############################
###Preliminaries
#set working directory
setwd("D:/Projects-1/4. Date compilation/Date calibration")
#removing old files
rm(list=ls())
#set libraries
library("Bchron")
library("rlist")
library("Deducer")
library("zoo")
library("SiZer")
library("geosphere")
library("Hmisc")
library("maptools")
library("raster")
library("sp")
library("rgdal")
library("gstat")
library("ape")
###########################################################################################################
#Function installation (check current)
###########################################################################################################
###Function for re-sampling date distributions. Data is a BChronCalibrate object, nsamp is number of re-samples.
SampleDates<-function(data,nsamp){
  matrix<-matrix(nrow=length(data),ncol=nsamp)
  for (i in 1:nsamp) {
    res<-as.integer(unlist(
      lapply(
        data,function(x){
          ages<-x$ageGrid
          densities<-x$densities
          samples<-sample(ages,1,prob=densities,replace=TRUE)
        }
      ),use.names=F))
    matrix[,i]<-res
  }
  return (matrix)
}
###########################################################################################################
###Function for re-sampling and comparing two sets of dates using a permutation t-test, nperm is number permutations for t-test.
CompareSampleDates<-function(data1,data2,nsamp,nperm){
  matrix1<-SampleDates(data1,nsamp)
  matrix2<-SampleDates(data2,nsamp)
  t.values<-rep(0,nsamp)
  p.values<-rep(0,nsamp)
  for (i in 1:nsamp){
    t.test<-perm.t.test(matrix1[,i],matrix2[,i],B=nperm)
    p.values[i]<-t.test$p.value
    t.values[i]<-t.test$statistic
  }
  print(noquote(c("mean p value=",round((mean(p.values)),3))))
  print(noquote(c("mean t value=",round((mean(t.values)),3))))
  output<-data.frame(t.values,p.values)
  return(output)
}
###########################################################################################################


#Reference. All dates

###Importing data
#Dates
Alldates<-read.csv("Alldates.csv")
#Associating spatial data
coordinates(Alldates)<-cbind(Alldates$E,Alldates$N)
proj4string(Alldates)<-CRS('+proj=longlat +datum=WGS84')

#Calibrating
Cal_dates_All<-BchronCalibrate(Alldates$Basal_date,Alldates$Error,Alldates$Cal_curve,Alldates$Site_name)
###Sampling dates, 1000 samples
Alldates_sample<-SampleDates(Cal_dates_All,1000)
###Mean and output
Alldates_average<-setNames(c("Alldates",(mean(colMeans(Alldates_sample))),(quantile(colMeans(Alldates_sample), c(0.05,0.95)))),c("approach","mean","5%","95%"))


###########################################################################################################
#APPROACH 1. Multiple dates per site
#Options: averaging by site, only dates where multiple dates per site. 

###Importing data
#Dates
multidates<-read.csv("multidates.csv")
#Sites as factor
multidates$Site<-as.factor(multidates$Site)

###Calibrating dates using Bchron
Cal_dates_multi<-BchronCalibrate(multidates$Basal_date,multidates$Error,multidates$Cal_curve,ids=multidates$Site_name)
##SPD<-need to load function
#Cal_dates_multiSPD<-SPD(Cal_dates_multi)
#plot(Cal_dates_multiSPD)

###Sampling dates, 1000 samples
multidates_sample<-SampleDates(Cal_dates_multi,1000)

###All dates from multi-datedsites
colmeans_multidates_sample<-colMeans(multidates_sample)#average age per rep
rowmeans_multidates_sample<-rowMeans(multidates_sample)#average for each dated point
hist(rowmeans_multidates_sample)
multidate_average<-setNames(c("multidateaverage",(mean(colmeans_multidates_sample)),(quantile(colmeans_multidates_sample, c(0.05,0.95)))),c("approach","mean","5%","95%"))

###By sites: mean, min, max
Sitemeans_multidates_sample<-apply(multidates_sample, 2, function (x) tapply(x, multidates$Site, mean))
Sitemeans_multidates_sample_rowMeans<-rowMeans(Sitemeans_multidates_sample)#average mean for each site
hist(Sitemeans_multidates_sample_rowMeans)
Sitemeans_average<-setNames(c("Bysitemeans",(mean(colMeans(Sitemeans_multidates_sample))),(quantile(colMeans(Sitemeans_multidates_sample), c(0.05,0.95)))),c("approach","mean","5%","95%"))
rm(Sitemeans_multidates_sample)#removing full file to save space
Sitemins_multidates_sample<-apply(multidates_sample, 2, function (x) tapply(x, multidates$Site, min))
Sitemins_multidates_sample_rowMeans<-rowMeans(Sitemins_multidates_sample)#average min for each site
hist(Sitemins_multidates_sample_rowMeans)
Sitemins_average<-setNames(c("Bysitemins",(mean(colMeans(Sitemins_multidates_sample))),(quantile(colMeans(Sitemins_multidates_sample), c(0.05,0.95)))),c("approach","mean","5%","95%"))
rm(Sitemins_multidates_sample)#removing full file to save space
Sitemaxs_multidates_sample<-apply(multidates_sample, 2, function (x) tapply(x, multidates$Site, max))
Sitemaxs_multidates_sample_rowMeans<-rowMeans(Sitemaxs_multidates_sample)#average max for each site
hist(Sitemaxs_multidates_sample_rowMeans)
Sitemaxs_average<-setNames(c("Bysitemaxs",(mean(colMeans(Sitemaxs_multidates_sample))),(quantile(colMeans(Sitemaxs_multidates_sample), c(0.05,0.95)))),c("approach","mean","5%","95%"))
rm(Sitemaxs_multidates_sample)#removing full file to save space
rm(multidates_sample)#removing full file to save space


###########################################################################################################
#APPROACH 2. Per grid-cell: oldest, mean, youngest
#NB in Holmquist this is 2808 dates for 1017 cells, c2.7 dates per cell. Here 267 dates for 73 non-blank cells. 

##########################################################################################################
#FUNCTION DEFINITION: RasterDates

RasterDates<-function(data,nsamp=1000,easting,northing,name,raster.cols=10,raster.rows=25,func=mean,raster.cutoff=4){
  raster_mean_output<-rep(0,nsamp)
  raster_results<-matrix(nrow=raster.cols*raster.rows,ncol=nsamp)
  for (i in 1:1000){
    Working<-data.frame(easting,northing,data[,i])
    coordinates(Working)<-cbind(Working[,1],Working[,2])
    proj4string(Working)<-CRS('+proj=longlat +datum=WGS84')
    blankraster<-raster(extent(Working), ncols=raster.cols, nrows=raster.rows)
    raster<-rasterize(Working, blankraster, data[,i], fun=function(x,na.rm){ifelse(length(x)<4,NA,func(x))})
    raster_mean<-mean(getValues(raster),na.rm=T)
    raster_mean_output[i]<-raster_mean
    raster_results[,i]<-getValues(raster)
  }
  raster_results_rowMeans<-rowMeans(raster_results)#this is the mean for each grid cell
  hist(raster_results_rowMeans)
  Raster_average_summary<-setNames(c(name,(mean(raster_mean_output)),(quantile(raster_mean_output, c(0.05,0.95)))),c("approach","mean","5%","95%"))
  output<-list("results"=raster_results_rowMeans,"summary"=Raster_average_summary)
  return(output)
}

###Per grid cell runs.
#Max
RasterMax<-RasterDates(data=Alldates_sample,nsamp=1000,easting=Alldates$E,northing=Alldates$N,name="Rastermaxall",raster.cols=10,raster.rows=25,func=max,raster.cutoff=4)
#Min
RasterMin<-RasterDates(data=Alldates_sample,nsamp=1000,easting=Alldates$E,northing=Alldates$N,name="Rasterminall",raster.cols=10,raster.rows=25,func=min,raster.cutoff=4)
#Mean
RasterMean<-RasterDates(data=Alldates_sample,nsamp=1000,easting=Alldates$E,northing=Alldates$N,name="Rastermeannall",raster.cols=10,raster.rows=25,func=mean,raster.cutoff=4)

###########################################################################################################
#APPROACH 3. Kriging
#

#Moran's I
DistanceMatrix<-as.matrix(dist(cbind(Alldates$E,Alldates$N)))  #assuming flat earth
DistanceMatrixInv<-1/DistanceMatrix
diag(DistanceMatrixInv)<-0
Moran.I(Alldates$Basal_date,DistanceMatrixInv)


#variogram
date.variogram<-variogram(Alldates$Basal_date~1, location=coordinates(Alldates), Alldates)
plot(date.variogram)
variogram.fit<-fit.variogram(date.variogram,model=vgm())#unfinished







###########################################################################################################
#APPROACH 4. Stacking following holmquist
#Seemingly not relevant because involves sub-setting the globe by regions and then scaling peatland area accordingly. Not really relevant to the UK?
#I guess? Or worth trying on a regional level?


#Calculate seperate spds for England scotland wales, weight by area, take mean...?

################
#Testing for differences- modify comparesampledates script to do this. Compare one 'sample dates' object with one other series derived from approaches above.
#for this initial test using Alldates but to be replaced with actual basal dates data
#1. comparing multidates with X (previous function)
MultidateComparison<-CompareSampleDates(Cal_dates_All,Cal_dates_multi,1000,1000)#replace with correct file

#2. By sites. Modified function.
CompareDatesNoSamp<-function(data1,data2,nsamp,nperm){
  matrix1<-data1
  matrix2<-SampleDates(data2,nsamp)
  t.values<-rep(0,nsamp)
  p.values<-rep(0,nsamp)
  for (i in 1:nsamp){
    t.test<-perm.t.test(matrix1,matrix2[,i],B=nperm)
    p.values[i]<-t.test$p.value
    t.values[i]<-t.test$statistic
  }
  print(noquote(c("mean p value=",round((mean(p.values)),3))))
  print(noquote(c("mean t value=",round((mean(t.values)),3))))
  output<-data.frame(t.values,p.values)
  return(output)
}
  #Applying
SitemeanCompare<-CompareDatesNoSamp(Sitemeans_multidates_sample_rowMeans,Cal_dates_All,1000,1000)
SiteminCompare<-CompareDatesNoSamp(Sitemins_multidates_sample_rowMeans,Cal_dates_All,1000,1000)
SitemaxCompare<-CompareDatesNoSamp(Sitemaxs_multidates_sample_rowMeans,Cal_dates_All,1000,1000)

#3. Rastered.
RasterMaxCompare<-CompareDatesNoSamp(RasterMax$results[!is.na(RasterMax$results)],Cal_dates_All,1000,1000)
RasterMinCompare<-CompareDatesNoSamp(RasterMin$results[!is.na(RasterMin$results)],Cal_dates_All,1000,1000)
RasterMeanCompare<-CompareDatesNoSamp(RasterMean$results[!is.na(RasterMean$results)],Cal_dates_All,1000,1000)



#################################################################################################################
###Summary file
SiteMethodsSummary<-rbind(Alldates_average,multidate_average,Sitemeans_average,Sitemins_average,Sitemaxs_average, RasterMean$summary, RasterMax$summary, RasterMin$summary)
rm(Alldates_average,multidate_average,Sitemeans_average,Sitemins_average,Sitemaxs_average)
