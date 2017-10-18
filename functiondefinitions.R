###Function definitions for RC project
###############################
###Preliminaries
#set working directory
setwd("D:/Projects-1/4. Date compilation/Date calibration")
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
###Create a summed probability density curve from a BChronCalibrate dataset
SPD<-function(data){
            Dates<- -2000:70000
            SummedDensity<-rep(0, length(Dates))
            for (i in 1:length(data)) {
              Matches<-match(data[[i]]$ageGrid, Dates)
              SummedDensity[Matches]=SummedDensity[Matches]+data[[i]]$densities
              }
            Output<-data.frame(Dates,SummedDensity)
            Output<-Output[min(which(Output$SummedDensity!=0)):max(which(Output$SummedDensity!=0)),]
            return(Output)
            }


###########################################################################################################
###'Uncalibrate' function from Crema et al. (2016 PLoS ONE 11(4): e0154809)
uncalibrate<-function(dates,error,calCurves='intcal13',random=TRUE){
                require(Bchron) # Bchron v4.0 
                pathToCalCurves=system.file("data", package = "Bchron")
                calCurveFile = paste(pathToCalCurves, "/", calCurves,".txt.gz", sep = "")
                calcurve=as.matrix(read.table(calCurveFile))[,1:3]
                colnames(calcurve) <- c("CALBP", "C14BP", "Error")
  
                ## uncalibrate CAL BP dates, interpolating with approx
                dates <- data.frame(approx(calcurve, xout = dates))
                colnames(dates) <- c("CALBP", "C14BP")
                calcurve.error <- approx(calcurve[,c(1,3)], xout = dates$CALBP)$y
                dates$Error <- sqrt(error^2 + calcurve.error^2)
                if(random==TRUE){dates$C14.Age=round(rnorm(nrow(dates),mean=dates$C14BP,sd=dates$Error))}
                return(dates)
                }

###########################################################################################################
###Generate random summed probability distributions
###Generates random calendar dates evenly distributed across a range selected from a 'SampleDates' matrix, samples errors from an observed distribution, uncalibrates using code from Crema et al. (2016) and re-calibrates using BchronCalibrate, outputting a summed probability density

GenRand<-function(ndates,error,limits){
          min<-min(limits[,sample(ncol(limits), 1)])
          max<-max(limits[,sample(ncol(limits), 1)])
          randcaldates<-sample(min:max, ndates, replace=T)
          randerrors<-sample(error, ndates, replace=T)
          uncalranddates<-uncalibrate(randcaldates,randerrors)
          curve<-rep("intcal13",ndates)
          recalranddates<-BchronCalibrate(uncalranddates$C14BP,uncalranddates$Error,curve)
          RandSPD<-SPD(recalranddates)
          return(RandSPD)
          }

###NB check how Crema treat errors

###########################################################################################################
###Generates confidence intervals from randomly generated SPDs
###nsamps is number of random iterations

RandCI<-function(ndates,error,limits,nsamps){
          Dates<- -2000:70000
          matrix<-matrix(nrow=length(Dates),ncol=nsamps)
          matrix[]<-0L
          rownames(matrix)<-Dates
          for (i in 1:nsamps){
             randrep<-GenRand(ndates=ndates,error=error,limits=limits)
             Matches<-match(randrep$Dates, Dates)
             matrix[,i][Matches]<-matrix[,i][Matches]+randrep$SummedDensity
             }
          lowerCI<-apply(matrix,1,quantile,prob=0.025)
          upperCI<-apply(matrix,1,quantile,prob=0.975)        
          CIs<-data.frame(Dates,lowerCI,upperCI)
          CIs<-CIs[min(which(CIs$lowerCI!=0)):max(which(CIs$lowerCI!=0)),]
          return(CIs)
          } #Need to think about smoothing

###########################################################################################################
###Uses a re-sampling approach to test for difference between SPDs. 
###Method correlates two SPDs and then compares correlation coefficient (Spearman Rs as default) to correlation with random shuffling of dates between two datasets of equivalent size to original
###Parameters: names/dates/errors/calcurves for two datasets, nperm=number cycles, cormethod=correlation coefficient: "spearman" or "pearson"

SPDDiff<-function(nperm, dates1, error1, calcurve1, names1, dates2, error2, calcurve2, names2, cormethod="spearman"){
          Cal1<-BchronCalibrate(dates1, error1, calcurve1, names1)
          Cal2<-BchronCalibrate(dates2, error2, calcurve2, names2)
          SPD1<-SPD(Cal1)
          SPD2<-SPD(Cal2) 
          #Matching, combining and correlating the two SPDs
          Dates<- -2000:70000
          matrix<-matrix(nrow=length(Dates),ncol=2)
          matrix[]<-0L
          rownames(matrix)<-Dates
          colnames(matrix)<-c("SPD1","SPD2")
          Matches1<-match(SPD1$Dates, Dates)
          Matches2<-match(SPD2$Dates, Dates)
          matrix[,1][Matches1]<-matrix[,1][Matches1]+SPD1$SummedDensity
          matrix[,2][Matches2]<-matrix[,2][Matches2]+SPD2$SummedDensity
          zerocheck<-matrix[,1]+matrix[,2]
          matrix<-matrix[min(which(zerocheck!=0)):max(which(zerocheck!=0)),]
          Cor<-rcorr(matrix, type=cormethod)
          #Establishing an output vector and starting a loop
          CorrelationResults<-rep(0,nperm)
          for (i in 1:nperm){
              #Combine to single matrix
              set1<-data.frame(names1, dates1, error1, calcurve1, stringsAsFactors=F)
              colnames(set1)<-c("names","dates","error","calcurve")
              set2<-data.frame(names2, dates2, error2, calcurve2, stringsAsFactors=F)
              colnames(set2)<-c("names","dates","error","calcurve")
              combinedDF<-rbind(set1,set2)
              #Re-sample into sets of same size as original series
              ndates<-length(dates1)
              indices<-sample(nrow(combinedDF), ndates, replace=F)
              randsamp<-combinedDF[indices, ]
              randrest<-combinedDF[-indices, ]
              #Re-calibrate and calculate SPDs
              Calrand<-BchronCalibrate(randsamp$dates, randsamp$error, randsamp$calcurve, randsamp$names)
              Calrest<-BchronCalibrate(randrest$dates, randrest$error, randrest$calcurve, randrest$names)
              randSPD<-SPD(Calrand)
              restSPD<-SPD(Calrest)
              #Combine and then compare SPDs based on randomised data
              RandDates<- -2000:70000
              Randmatrix<-matrix(nrow=length(RandDates),ncol=2)
              Randmatrix[]<-0L
              rownames(Randmatrix)<-RandDates
              colnames(Randmatrix)<-c("SPD1","SPD2")
              RandMatches1<-match(randSPD$Dates, RandDates)
              RandMatches2<-match(restSPD$Dates, RandDates)
              Randmatrix[,1][RandMatches1]<-Randmatrix[,1][RandMatches1]+randSPD$SummedDensity
              Randmatrix[,2][RandMatches2]<-Randmatrix[,2][RandMatches2]+restSPD$SummedDensity
              Randzerocheck<-Randmatrix[,1]+Randmatrix[,2]
              Randmatrix<-Randmatrix[min(which(Randzerocheck!=0)):max(which(Randzerocheck!=0)),]
              RandCor<-rcorr(Randmatrix, type=cormethod)
              CorrelationResults[i]<-RandCor$r[1,2]
              }
          P<-(length(CorrelationResults[CorrelationResults>Cor$r[1,2]]))/nperm
          print(noquote(c("P=",1-P)))
          output<-list(CorrelationCoefficient=Cor$r[1,2],RandomCorrelations=CorrelationResults)
          return(output)
          }

###########################################################################################################
###Uses bootstrapping to derive confidence intervals for an SPD. 
###Method selects n samples with replacements, repeats multiple times and calculates 95% confidence intervals based on the results. 
###Parameters: names/dates/errors/calcurves for the dataset, nperm=number bootstrap cycles


SPDboot<-function(nboot, dates, error, calcurve, names){
          data<-data.frame(names, dates, error, calcurve, stringsAsFactors=F)
          colnames(data)<-c("names","dates","error","calcurve")
          Dates<- -2000:70000
          matrix<-matrix(nrow=length(Dates),ncol=nboot)
          matrix[]<-0L
          rownames(matrix)<-Dates
          for (i in 1:nboot){
              selected<-sample(nrow(data), nrow(data), replace=T)
              sample<-data[selected, ]
              Calibrated<-BchronCalibrate(sample$dates, sample$error, sample$calcurve, sample$names)
              SPD<-SPD(Calibrated)
              matches<-match(SPD$Dates, Dates)
              matrix[,i][matches]<-matrix[,i][matches]+SPD$SummedDensity
              }
          lowerCI<-apply(matrix,1,quantile,prob=0.025)
          upperCI<-apply(matrix,1,quantile,prob=0.975)        
          CIs<-data.frame(Dates,lowerCI,upperCI)
          CIs<-CIs[min(which(CIs$lowerCI!=0)):max(which(CIs$lowerCI!=0)),]
          return(CIs)
          }## double check all this when fresh
          ### also need to think about smoothing


###########################################################################################################
###Transfers radiocarbon dates to raster cells, using re-sampling of date distribution to encompass non-linearity.
###Method takes an output from SampleDates, assigns coordinates, fits to specified raster grid and calculates values using specified function with specified cut-off. Assumes WGS84 coordinates. Returns a summary file and output based on date means. 
###Parameters: data is a SampleDates output with assocated coordinates, raster specification, function to derive values and cut-off for minimum number of points informing a grid cell. 


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



###############################################################################################################
#### Note that although the code in this file is functional (to the best of my knowledge) it is is not      ###
#### neccessarily efficient in terms of time or memory usage. Some functions are very slow to run and       ###
#### many could probably be improved with better coding. Use at your own risk.                              ###
###############################################################################################################
