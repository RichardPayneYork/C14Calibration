###Initial script for RC project
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
###Importing data
#All dates
Alldates<-read.csv("alldates.csv")
#'A' dates
Adates<-read.csv("adates.csv")
###Calibrating dates using Bchron
Cal_dates<-BchronCalibrate(Adates$Basal_date,Adates$Error,Adates$Cal_curve,ids=Adates$Site_name)
Cal_datesAll<-BchronCalibrate(Alldates$Basal_date,Alldates$Error,Alldates$Cal_curve,ids=Alldates$Site_name)
###Comparing calibrated dates using re-sampling approach
#Data set 1
#Set number of permutations and file to use
nperm=3
dat=Cal_datesAll
#Creating a data matrix, re-sampling date probability distributions and outputting to matrix
count<-length(dat)
matrix<-matrix(nrow=count,ncol=nperm)
for (i in 1:nperm) {
  res<-as.integer(unlist(
    lapply(
      dat,function(x){
      ages<-x$ageGrid
      densities<-x$densities
      samples<-sample(ages,1,prob=densities,replace=TRUE)
      }
    ),use.names=F))
  matrix[,i]<-res
}
###Comparing calibrated dates using re-sampling approach
#Data set 2
#Set number of permutations and file to use
nperm=1000
dat2=Cal_dates
#Creating a data matrix, re-sampling date probability distributions and outputting to matrix
count2<-length(dat2)
matrix2<-matrix(nrow=count2,ncol=nperm)
for (i in 1:nperm) {
  res<-as.integer(unlist(
    lapply(
      dat2,function(x){
        ages<-x$ageGrid
        densities<-x$densities
        samples<-sample(ages,1,prob=densities,replace=TRUE)
      }
    ),use.names=F))
  matrix2[,i]<-res
}
####Using permutation t-test to test for difference between two matrices above and outputting to vectors
t.values<-rep(0,nperm)
p.values<-rep(0,nperm)
for (i in 1:nperm){
  t.test<-perm.t.test(matrix[,i],matrix2[,i],B=nperm)
  print(i)
  p.values[i]<-t.test$p.value
  t.values[i]<-t.test$statistic
}
mean(t.values)
mean(p.values)
boxplot(t.values)
boxplot(p.values)
####Deleting these large matrices when finished
rm(matrix)
rm(matrix2)
####Using Bchron to calculate summed density
#All data
DensityAll<-BchronDensity(Alldates$Basal_date,Alldates$Error,Alldates$Cal_curve)
#####These lines simply tweak and apply the existing plotting functionality in Bchron so that the summed 
#####density results are output to a vector for off-line plotting. 
plot.Bchron<-function (x, plotDates = TRUE, plotSum = FALSE, ...) 
{
  n = length(x$calAges)
  thetaRange = range(x$calAges[[1]]$ageGrid)
  for (i in 2:n) thetaRange = range(c(thetaRange, x$calAges[[i]]$ageGrid))
  dateGrid = seq(round(thetaRange[1] * 0.9, 3), round(thetaRange[2] * 
                                                        1.1, 3), length = 1000)
  gauss <- function(x, mu, sig) {
    u <- (x - mu)/sig
    y <- exp(-u * u/2)
    y
  }
  gbase <- function(x, mus) {
    sig <- (mus[2] - mus[1])/2
    G <- outer(x, mus, gauss, sig)
    G
  }
  Gstar = gbase(dateGrid, x$mu)
  dens = vector(length = length(dateGrid))
  for (i in 1:nrow(x$p)) {
    dens = dens + Gstar %*% x$p[i, ]
  }
  densFinal = dens/sum(dens)
  graphics::plot(dateGrid, densFinal, type = "l", ylab = "Density", 
                 ylim = range(c(0, densFinal)), ...)
  if (plotDates) {
    yHeight = graphics::par("usr")[4]
    myCol = grDevices::rgb(190/255, 190/255, 190/255, 0.4)
    for (i in 1:n) {
      graphics::polygon(x$calAges[[i]]$ageGrid, 0.3 * yHeight * 
                          x$calAges[[i]]$densities/max(x$calAges[[i]]$densities), 
                        col = myCol, border = NA)
    }
  }
  if (plotSum) {
    yHeight = graphics::par("usr")[4]
    thetaRange = range(x$calAges[[1]]$ageGrid)
    for (i in 2:n) thetaRange = range(c(thetaRange, x$calAges[[i]]$ageGrid))
    dateGrid = seq(round(thetaRange[1] * 0.9, 0), round(thetaRange[2] * 
                                                          1.1, 0), by = 1)
    sumDens = rep(0, length(dateGrid))
    for (i in 1:n) {
      matchRows = match(x$calAges[[i]]$ageGrid, dateGrid)
      sumDens[matchRows] = sumDens[matchRows] + x$calAges[[i]]$densities
      if (any(is.na(matchRows))) 
                stop()
         }
     graphics::lines(dateGrid, sumDens * yHeight/max(sumDens), 
                    col = "red")
  }
  dens<-(sumDens)
  age<-(dateGrid)
  data<-cbind(age,dens)
}
DensityPlotAll<-plot.Bchron(DensityAll,plotSum=TRUE)
plot##Applying 500 year moving mean smoother to this summed density plot
plot(DensityPlotAll[,1],DensityPlotAll[,2])
smoothed<-rollapply(DensityPlotAll,500,mean)
plot(smoothed)
##Alternatively: using smoothing spline to smooth this density plot, spar value determines smoothing
smoother<-smooth.spline(DensityPlotAll[,1],DensityPlotAll[,2],spar=1)
plot(smoother)
###Applying SiZer to density plot. Based on grid length of 500 intervals and first derivative. This is slow. 
datasizerbig5<-SiZer(DensityPlotAll[,1],DensityPlotAll[,2],h=c(.5,2000),x.grid=1000,derv=1)
plot(datasizerbig5, colorlist = c("red", "purple", "blue", "black"))
?SiZer
str(datasizerbig5)
typeof(datasizerbig5)
write.table(datasizerbig5$x.grid,file="sizergrid.txt")
write.table(datasizerbig5$h.grid,file="sizerhgrid.txt")
write.table(datasizerbig5$slopes,file="sizerslopes.txt")
plot.SiZer
####re-writing sizer plot function to remove white line
sizerplot<-function (x, ylab = expression(log[10](h)), colorlist = c("red", 
                                                          "purple", "blue", "grey"), ...) 
{
  temp <- factor(x$slopes)
  final.colorlist <- NULL
  if (is.element("-1", levels(temp))) 
    final.colorlist <- c(final.colorlist, colorlist[1])
  if (is.element("0", levels(temp))) 
    final.colorlist <- c(final.colorlist, colorlist[2])
  if (is.element("1", levels(temp))) 
    final.colorlist <- c(final.colorlist, colorlist[3])
  if (is.element("2", levels(temp))) 
    final.colorlist <- c(final.colorlist, colorlist[4])
  temp <- matrix(as.integer(factor(x$slopes)), nrow = dim(x$slopes)[1])
  image(x$x.grid, log(x$h.grid, 10), t(temp), col = final.colorlist, 
        ylab = ylab, ...)
  x.midpoint <- diff(range(x$x.grid))/2 + min(x$x.grid)
  }
###applying re-written sizer plot function
sizerplotbig5better<-sizerplot(datasizerbig5)
####
Cal_dates$ageGrid
sum(Cal_dates == 7000)
un<-unlist(Cal_dates)
list.count(Cal_dates$Dod$ageGrid,4000)
typeof(un)
log10(2000)
#rm(list = ls(all.names = TRUE))

##############################################################################################
#############Using archaeo dates to look at similarity/difference to basal peat dates
#'Archaeo' dates
Archdates<-read.csv("ArchaeoDates.csv")
###Calibrating dates using Bchron
Cal_datesArch<-BchronCalibrate(Archdates$Date,Archdates$Error,Archdates$Cal_curve,ids=Archdates$ID)
###Using Bchron to calculate summed density- this is SLOW!
###This hits against memory limits therefore need to increase before running
memory.limit(size=20000)
DensityArch<-BchronDensity(Archdates$Date,Archdates$Error,Archdates$Cal_curve)
###Applying tweaked BChron density plotting function
DensityPlotArch<-plot.Bchron(DensityArch,plotSum=TRUE)
plot(DensityPlotArch[,1],DensityPlotArch[,2],xlim=c(0,10000))
###Comparing calibrated dates using re-sampling approach
#Data set 3: Archaeological RC dates
#Set number of permutations and file to use
nperm=1000
dat3=Cal_datesArch
#Creating a data matrix, re-sampling date probability distributions and outputting to matrix
count<-length(dat3)
matrix3<-matrix(nrow=count,ncol=nperm)
for (i in 1:nperm) {
  res<-as.integer(unlist(
    lapply(
      dat3,function(x){
        ages<-x$ageGrid
        densities<-x$densities
        samples<-sample(ages,1,prob=densities,replace=TRUE)
      }
    ),use.names=F))
  matrix3[,i]<-res
}
##Using permutation t-test to test for difference between two matrices above and outputting to vectors
t.values<-rep(0,nperm)
p.values<-rep(0,nperm)
for (i in 1:nperm){
  t.test<-perm.t.test(matrix[,i],matrix3[,i],B=nperm)
  print(i)
  p.values[i]<-t.test$p.value
  t.values[i]<-t.test$statistic
}
mean(t.values)
mean(p.values)
boxplot(t.values[1:1000])
boxplot(p.values [1:1000])
mean(matrix3)
mean(matrix)
###p-values all zero. maybe due to inclusion of pre-holocene dates? consider subsetting?
#############################################################################################################
##Calculating spatio-temporal distance between peat dates and nearest 
#Spatial distance with Vincenty Ellipsoid method [slow!]
Arch<-data.frame(Archdates$Long,Archdates$Lat)
Peat<-data.frame(Alldates$E,Alldates$N)
dist<-distm(Arch[,c('Archdates.Long','Archdates.Lat')], Peat[,c('Alldates.E','Alldates.N')], fun=distVincentyEllipsoid)
#Centering and scaling. 
dist=(dist-mean(dist))/(sd(dist))
###Temporal distance
#For first column only calculates a dissimilarity matrix
count1<-nrow(matrix)
count2<-nrow(matrix3)
output<-array(dim=c(count1,count2))
for (j in 1:count2){
for (i in 1:count1){
  output[i,j]<-sqrt((matrix[i,1]-matrix3[j,1])^2)
    }
}
#Next stage- z-scores, multiply up by spatial distance then loop for other columns 
#centering, scaling and transposing
output=t((output-mean(output))/(sd(output)))
#making values positive and multiplying spatial and using pythagoras theorem to give spatio-temporal distance
#fulldist=sqrt(((abs(dist))^2)+((abs(output))^2))
fulldist=sqrt(((abs(dist)))+((abs(output))))
###finding cell with smallest combined distance
mins<-rep(0,count1)
for (k in 1:count1){
  mins[k]<-which.min(fulldist[,k])
}
mins 
###returning absolute values
nearestdates<-matrix3[mins,1]
##need to double check this
results<-data.frame(nearestdates,matrix[,1])
plot(results)
boxplot(nearestdates)
boxplot(matrix[,1])
rcorr(nearestdates,matrix[,1])


str(matrix[,1])
str(nearestdates)

Alldates[266,]
Archdates[1179,]
#something wrong here- giving odd-looking results. Seems to work for geog and temp distance seperately but not combined. More digging needed.


###Idea for comparing summed probability distributions between different sets of dates. 
###1. correlate both series using eg spearman
###2. combine both datasets and re-sample to give two new datasets of equivalent size (or re-sample with replacement from larger dataset?)
###3. sum probabilities, test correlation for each, repeat lots of times
###4. compare real correlation to simulations
###sensible???
###following Shennan et al 2013 could do cross correlation and take the single largest correlation coefficient as index

###Deviations in summed probability curve.
###Assume a null distribution of no change over time. Sample calendar dates randomly from this distribution.
###'Uncalibrate' using function in 'src.R' (need to check this), assign an error (pick randomly from real distribution?) and recalibrate. Repeat.
###Calculate 95% confidence intervals of these distributions and compare to real distribution to highlight excursions. May need to calculate z scores. Alternatively could just compare to *all* random simulations.  
# randomSDs<-sample(size=length(randomDates),error,replace=TRUE)
# simDates<-round(uncalibrate(randomDates,randomSDs,random=TRUE)[,2:3])

#importing intcal13 curve
calCurveFile<-paste("C:/Users/Richard/Documents/R/win-library/3.1/Bchron/data/intcal13.txt.gz")
calcurve<-as.matrix(read.table(calCurveFile))[,1:3]
colnames(calcurve)<- c("CALBP", "C14BP", "Error")

## uncalibrate CAL BP dates, interpolating with approx
dates <- data.frame(approx(calcurve, xout = dates))
colnames(dates) <- c("CALBP", "C14BP")
calcurve.error <- approx(calcurve[,c(1,3)], xout = dates$CALBP)$y
dates$Error <- sqrt(error^2 + calcurve.error^2)
if(random==TRUE){dates$C14.Age=round(rnorm(nrow(dates),mean=dates$C14BP,sd=dates$Error))}
return(dates)

#uncalibrate function from Crema et al. 2016
uncalibrate<-function(dates,error,calCurves='intcal13',random=TRUE)
{
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

###Tweaking re-sampling function to output min and max
CSamp<-SampleDates(Cal_dates,nsamp=100)


#Generate 267 random numbers within range of randomised distribution (may need to adjust to only use part of range)
##re-working
#'limits'uses the output of a SampleDates run to specify range of simulated dates based on a randomly selected run


GenRand<-function(nsamp,error,limits){
min<-min(limits[,sample(ncol(limits), 1)])
max<-max(limits[,sample(ncol(limits), 1)])
randcaldates<-sample(min:max, nsamp, replace=T)
randerrors<-sample(error, nsamp, replace=T)
uncalranddates<-uncalibrate(randcaldates,randerrors)
curve<-rep("intcal13",nsamp)
recalranddates<-BchronCalibrate(uncalranddates$C14BP,uncalranddates$Error,curve)
RandSPD<-SPD(recalranddates)
return(RandSPD)
}

Alldates$Cal_curve

XXX<-GenRand(nsamp=100,error=Alldates$Error,limits=CSamp)
plot(XXX)

error=Alldates$Error
limits=CSamp
nsamp=999
limits=NULL
nsamp=NULL
error=NULL
typeof(curve)
curve<-as.vector(curve)

randcalendardates<-sample(min(matrix):max(matrix), count1, replace=T)
randerrors<-sample(Alldates$Error, count1, replace=T)
uncalibratedrandomdates<-uncalibrate(randcalendardates,randerrors)





randdensity<-BchronDensityFast(uncalibratedrandomdates$C14BP,uncalibratedrandomdates$Error,Alldates$Cal_curve)







randdensityplot<-plot.Bchron.Fast(randdensity, plotDates=F, plotSum=T)








str(randdensityplot)
typeof(randdensityplot)
plot(randdensityplot)
str(randdensity)
typeof(randdensity)


# function below


str(randdensity)
randdensity$calAges
##repeating modification to Bchron density for BchrondensityFast
plot.Bchron.Fast<-function (x, plotDates = TRUE, plotSum = FALSE, ...) 
{
  n = length(x$calAges)
  nclusters = length(x$clusterRange)
  chosen = which.max(x$out$BIC)
  clusterMeans = x$out$parameters$mean
  clusterSds = sqrt(x$out$parameters$variance$sigmasq)
  clusterProps = x$out$parameters$pro
  thetaRange = round(c(min(clusterMeans - 3 * clusterSds), 
                       max(clusterMeans + 3 * clusterSds)), 0)
  thetaSeq = thetaRange[1]:thetaRange[2]
  dens = matrix(NA, ncol = x$out$G, nrow = length(thetaSeq))
  for (i in 1:ncol(dens)) {
    dens[, i] = stats::dnorm(thetaSeq, mean = clusterMeans[i], 
                             sd = clusterSds[i])
  }
  graphics::plot(thetaSeq, dens %*% clusterProps, type = "l", 
                 ...)
  for (i in 1:ncol(dens)) graphics::lines(thetaSeq, clusterProps[i] * 
                                            dens[, i], lty = 2)
  if (plotDates) {
    yHeight = graphics::par("usr")[4]
    myCol = grDevices::rgb(190/255, 190/255, 190/255, 0.4)
    for (i in 1:n) {
      graphics::polygon(x$calAges[[i]]$ageGrid, 0.3 * yHeight * 
                          x$calAges[[i]]$densities/max(x$calAges[[i]]$densities), 
                        col = myCol, border = NA)
    }
  }
  if (plotSum) {
    yHeight = graphics::par("usr")[4]
    thetaRange = range(x$calAges[[1]]$ageGrid)
    for (i in 2:n) thetaRange = range(c(thetaRange, x$calAges[[i]]$ageGrid))
    dateGrid = seq(round(thetaRange[1] * 0.9, 0), round(thetaRange[2] * 
                                                          1.1, 0), by = 1)
    sumDens = rep(0, length(dateGrid))
    for (i in 1:n) {
      matchRows = match(x$calAges[[i]]$ageGrid, dateGrid)
      sumDens[matchRows] = sumDens[matchRows] + x$calAges[[i]]$densities
      if (any(is.na(matchRows))) 
        stop()
    }
    graphics::lines(dateGrid, sumDens * yHeight/max(sumDens), 
                    col = "red")
      }
  Sumdens<-(sumDens)
   age<-(dateGrid)
  return(dens)
  #data<-cbind(age,Sumdens)
  }
###looking at structure of bchrondensityfast

modbcdF<-function (ages, ageSds, calCurves, pathToCalCurves = system.file("data", 
                                                                 package = "Bchron"), dfs = rep(100, length(ages)), samples = 2000, 
          G = 30) 
{
  if (length(ages) != length(ageSds)) 
    stop("ages and 1-sigma errors must be same length")
  if (length(ages) != length(calCurves)) 
    stop("ages and Calibration curves must be same length")
  x = BchronCalibrate(ages = ages, ageSds = ageSds, calCurves = calCurves, 
                      pathToCalCurves = pathToCalCurves, dfs = rep(100, length(ages)))
  n = length(x)
  thetaBig = vector(length = n * samples)
  for (i in 1:n) thetaBig[((i - 1) * samples + 1):(i * samples)] = sample(x[[i]]$ageGrid, 
                                                                          size = samples, prob = x[[i]]$densities, replace = TRUE)
  mclustOutput = mclust::densityMclust(data = thetaBig, G = G)
  output = list(out = mclustOutput, calAges = x)
  class(output) = "BchronDensityRunFast"
  return(output)

}
############
play<-function(x){
dateGrid = 1:20000
sumDens = rep(0, length(dateGrid))
n = length(x$calAges)
for (i in 1:n) {
  matchRows = match(x$calAges[[i]]$ageGrid, dateGrid)
  sumDens[matchRows] = sumDens[matchRows] + x$calAges[[i]]$densities
}
}

###works
Dates<--70:20000
SummedDensity<-rep(0, length(Dates))
for (i in 1:length(Cal_datesAll)) {
  Matches<-match(Cal_datesAll[[i]]$ageGrid, Dates)
  SummedDensity[Matches] = SummedDensity[Matches] + Cal_datesAll[[i]]$densities
}

##re-writing as a function
#trying to make work for NAs- works with larger range
SPD<-function(x){
            Dates<--2000:60000
            SummedDensity<-rep(0, length(Dates))
            for (i in 1:length(x)) {
                Matches<-match(x[[i]]$ageGrid, Dates)
                SummedDensity[Matches]=SummedDensity[Matches]+x[[i]]$densities
                 }
            Output<-data.frame(Dates,SummedDensity)
            Output<-Output[min(which(Output$SummedDensity!=0)):max(which(Output$SummedDensity!=0)),]
            return(Output)
            }


CSPDA<-SPD(Cal_dates)
plot(CSPDA)

##############################################################################################################
####Note that although the code in this file is functional it is is not neccessarily efficient in terms of ###
####time or memory usage. Some functions are very slow to run and many could probably be improved with     ###
####better coding.                                                                                         ###
##############################################################################################################