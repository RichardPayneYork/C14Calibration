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





###############################################################################################################
#### Note that although the code in this file is functional (to the best of my knowledge) it is is not      ###
#### neccessarily efficient in terms of time or memory usage. Some functions are very slow to run and       ###
#### many could probably be improved with better coding. Use at your own risk.                              ###
###############################################################################################################
