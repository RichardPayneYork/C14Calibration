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
library("R.basic")
library("scales")


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
###Arguments: names/dates/errors/calcurves for two datasets, nperm=number cycles, cormethod=correlation coefficient: "spearman" or "pearson"

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
###Arguments: names/dates/errors/calcurves for the dataset, nperm=number bootstrap cycles


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
###Arguments: data is a SampleDates output with assocated coordinates, raster specification, function to derive values and cut-off for minimum number of points informing a grid cell. 


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

###################################################################################################################
###Function for batch-importing Bacon files (assumes all .txt files in working directory are Bacon output) and calculating binned mean accumulation, accumulation weighted by model precision and inferred carbon accumulation based on other dated sites. Returns a list of data files, csv files in working directory and a plot of C accumulation.
###Arguments: Max,min,bin specify bins used. Names of output files. Number of permutations used. Vectors of carbon density (carbon.data) and the age of these samples (carbon.age), if ommited calculates just depth accumulation. Whether C data selection is entirely random (by.depth=F) or by depth (by.depth=T). 

Accumulation<-function(carbon.age=NULL, carbon.data=NULL, nperm=1000, min=0, max=10500, bin.interval=500, by.depth=T, accumulation.name="AccumulationRateOutput", wt.accumulation.name="WeightedAccumulationRateOutput", c.accumulation.name="CarbonOutput"){
        ### Batch importing BACON output files. 
        file.list<-list.files(pattern= "*.txt")
        data.list<-lapply(file.list, FUN=read.table, header=TRUE)
        ## Accumulation over every bin
        bins<-seq(min,max,bin.interval)
        bins.mid<-seq((min+(bin.interval/2)),(max-(bin.interval/2)),bin.interval)
        bins.max<-bins[2:length(bins)]
        matrix<-matrix(nrow=length(bins),ncol=length(file.list))
        for (i in 1:length(file.list)) {
             interp.depth<- approx(data.list[[i]][,5],data.list[[i]][,1], "linear", xout=bins, rule=1)
             matrix[,i]<-interp.depth$y
         }
           accumulation<-diff(matrix,lag=1)#calculating differences
           accumulation.correct<-accumulation/bin.interval #accumulation cm/yr
           mean.accumulation<-rowMeans(accumulation.correct, na.rm = T)
           number.accumulation<-(length(file.list))-rowSums(is.na(accumulation))
           accu.output<-data.frame(bins.mid,mean.accumulation,number.accumulation)
           write.csv(accu.output, paste(accumulation.name,".csv"))
          ## Weighting this by precision of age-depth model.
          ## Calculate min-max range for each bin for each site. 
          lower<-matrix(nrow=length(bins.mid),ncol=length(file.list))
          for (i in 1:length(file.list)) {
              interp.depth<- approx(data.list[[i]][,5],data.list[[i]][,1], method = "linear", xout=bins.mid, rule = 1)
              interp.min<- approx(data.list[[i]][,1],data.list[[i]][,2], method = "linear", xout=interp.depth$y, rule = 1)
              lower[,i]<-interp.min$y
              }
          upper<-matrix(nrow=length(bins.mid),ncol=length(file.list))
          for (i in 1:length(file.list)) {
            interp.depth<- approx(data.list[[i]][,5],data.list[[i]][,1], method = "linear", xout=bins.mid, rule = 1)
            interp.max<- approx(data.list[[i]][,1],data.list[[i]][,3], method = "linear", xout=interp.depth$y, rule = 1)
            upper[,i]<-interp.max$y
            }
          range<-upper-lower#full range of a-d models
          ## Calculate Z scores by rows and rescale to 0:1
          range.z<-t(apply(range, 1, FUN=zscore, na.rm=TRUE))#z-scores by rows and transposed back
          range.z.inv<-t(apply(range.z, 1, FUN=rescale, to = c(1, 0)))#rescaled and inverted 1:0 (1 best, 0 worst) and transposed back
          weightings<-range.z.inv+0.01 #adding a small offset so that no models don't contribute at all
          weightings.total<-rowSums(weightings, na.rm=TRUE)#total of weightings variables
          weighted<-(accumulation*weightings) #weighting
          wt.mean.accumulation<-(rowSums(weighted, na.rm = T))/weightings.total #summing and correcting 
          wt.mean.accumulation.correct<-wt.mean.accumulation/bin.interval #correcting to annual
          wt.output<-data.frame(bins.mid,wt.mean.accumulation.correct,number.accumulation)
          write.csv(wt.output, paste(wt.accumulation.name,".csv"))
          if(missing(carbon.data)) {
            return(list(accumulation=accu.output,weighted.accumulation=wt.output)) #to calculate accumulation if carbon data not supplied
          } else {
          ## Assigning carbon values to give inferred accumulation
          rand.accum.output<-matrix(nrow=nrow(accumulation.correct),ncol=nperm)
          wt.rand.accum.output<-matrix(nrow=nrow(accumulation.correct),ncol=nperm)
          for (j in 1:nperm){
            if(by.depth==F | is.null(carbon.age)){  #here splits, first option if by.depth=F, random selection of C density values
            rand.c.data<-matrix(nrow=nrow(accumulation.correct),ncol=ncol(accumulation.correct))
          for (i in 1:nrow(accumulation)) {
            rand.c.data[i,]<-sample(carbon.data[!is.na(carbon.data)], length(rand.c.data[i,]), replace = TRUE)#NB this fills full matrix with non-NA randomly selected values
          }}else{         #alternatively, here selects based on age bins
            categories<-(findInterval(carbon.age,bins))*bin.interval #max of age bin
            categories.carbon<-cbind(categories,carbon.data) #a double vector of carbon density and age
            null.matrix<-round((matrix+bins-matrix)) #making a matrix with age bin in place of depth, excluding the zero row.
            null.matrix<-null.matrix[2:nrow(null.matrix),] #removes top row to only show bin max
          for (k in 1:length(bins.max)){ 
              length.samp<-length(null.matrix[which(null.matrix==bins.max[k])]) #number of valid entries in age bin
              matched.values<-categories.carbon[,2][categories.carbon[,1]==bins.max[k]] #available carbon density values in bin
              if(length(matched.values)>0){
                null.matrix[which(null.matrix==bins.max[k])]<-sample(matched.values,length.samp,replace=T) #replacing bin values with a randomly selected (with replacement) density value
              } else {
                null.matrix[which(null.matrix==bins.max[k])]<-NA #if no available values assigning NA
            }}  
          rand.c.data<-null.matrix
        }
        rand.accum<-(rand.c.data*accumulation.correct)*10000 
        rand.mean.accum<-rowMeans(rand.accum, na.rm=T)
        rand.accum.output[,j]<-rand.mean.accum
        wt.rand.accum<-(rand.c.data*wt.mean.accumulation.correct)*10000 
        wt.rand.mean.accum<-rowMeans(wt.rand.accum, na.rm=T)
        wt.rand.accum.output[,j]<-wt.rand.mean.accum
      }
      overall.rand.accum.output<-rowMeans(rand.accum.output, na.rm=T)#means
      wt.overall.rand.accum.output<-rowMeans(wt.rand.accum.output, na.rm=T)
      rand.c.output<-data.frame(bins.mid,overall.rand.accum.output,wt.overall.rand.accum.output,number.accumulation)
      write.csv(rand.c.output, paste(c.accumulation.name,".csv"))
      plot(bins.mid,overall.rand.accum.output)
      return(list(bins=bins.mid,accumulation=accu.output,weighted.accumulation=wt.output,carbon.accumulation=rand.c.output))
  }}


###########################################################################################################
###Function to calculate overall core accumulation rate. Assumes every text file in the working directory is a Bacon output. Only input is specified name of output file. 

SiteAccumulation<-function(site.accumulation.name="output"){
      file.list<-list.files(pattern = "*.txt")
      data.list<-lapply(file.list, FUN=read.table, header=TRUE)
      max.min<-matrix(nrow=length(file.list),ncol=2)
      for (i in 1:length(file.list)) {
        max.min[i,1]<-min(data.list[[i]][,5])
        max.min[i,2]<-max(data.list[[i]][,5])
      }
      age.span<-max.min[,2]-max.min[,1]
      max.depth<-rep(0,length(file.list))
      for (i in 1:length(file.list)) {
        max.depth[i]<-max(data.list[[i]][,1])
      }
      site.accum<-age.span/max.depth
      site.accum.output<-data.frame(file.list,site.accum)
      write.csv(site.accum.output, paste(site.accumulation.name,".csv"))
      return(site.accum.output)
}


###############################################################################################################
#### Note that although the code in this file is functional (to the best of my knowledge) it is is not      ###
#### neccessarily efficient in terms of time or memory usage. Some functions are very slow to run and       ###
#### many could probably be improved with better coding. Use at your own risk.                              ###
###############################################################################################################
