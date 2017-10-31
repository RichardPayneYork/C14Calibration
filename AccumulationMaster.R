##### New R script for project on FloCo carbon accumulation project.
# Emptying workspace
rm(list=ls(all=TRUE))
# Loading packages
library("R.basic")
library("scales")


CarbonFile<-read.csv("Carbondata.csv", header=T)
CarbonFile[CarbonFile=="NaN"]<-NA #replacing Joss's NaN with R's NA
CarbonFile<-na.omit(CarbonFile) #removing all the NAs
CarbonFile$carbon<-CarbonFile$carbon/100 #converting percentage C to proportion
CarbonFile$density<-CarbonFile$carbon*CarbonFile$bulk.density #calculating carbon density. carbon per cm2
Carbondata<-CarbonFile$density

###################################################################################################################
#Function for batch-importing Bacon files (assumes all .txt files in working directory are Bacon output) and calculating binned mean accumulation, accumulation weighted by model precision and inferred carbon accumulation based on other dated sites. Returns a list of data files, csv files in working directory and a plot of C accumulation.
#Arguments: Max,min,bin specify bins used. Names of output files. Number of permutations used. Vectors of carbon density (carbon.data) and the age of these samples (carbon.age). Whether C data selection is entirely random (by.depth=F) or by depth (by.depth=T). 

Accumulation<-function(carbon.age, carbon.data, nperm=1000, min=0, max=10500, bin.interval=500, by.depth=T, accumulation.name="AccumulationRateOutput", wt.accumulation.name="WeightedAccumulationRateOutput", c.accumulation.name="CarbonOutput"){
  ### Stage 1. Batch importing BACON output files. 
  file.list<-list.files(pattern= "*.txt")
  data.list<-lapply(file.list, FUN=read.table, header=TRUE)
  ## Stage 2. accumulation over every bin
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
  ## Stage 3. Weighting this by precision of age-depth model.
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
    ## Stage 4
    ## Assigning carbon values to give inferred accumulation
    rand.accum.output<-matrix(nrow=nrow(accumulation.correct),ncol=nperm)
    wt.rand.accum.output<-matrix(nrow=nrow(accumulation.correct),ncol=nperm)
    categories<-(findInterval(carbon.age,bins))*bin.interval #max of age bin
    categories.carbon<-cbind(categories,carbon.data) #a double vector of carbon density and age

    for (j in 1:nperm){
      if(by.depth==F | missing(carbon.age)){  #here splits, first option if by.depth=F, random selection of C density values
        rand.c.data<-matrix(nrow=nrow(accumulation.correct),ncol=ncol(accumulation.correct))
        for (i in 1:nrow(accumulation)) {
          rand.c.data[i,]<-sample(carbon.data[!is.na(carbon.data)], length(rand.c.data[i,]), replace = TRUE)#NB this fills full matrix with non-NA randomly selected values
        }}else{         #alternatively, here selects based on age bins
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
###########################################################################################################

output<-Accumulation(max=10500, min=0, bin.interval=500, accumulation.name="AccumulationRateOutput", wt.accumulation.name="WeightedAccumulationRateOutput", c.accumulation.name="CarbonOutput", by.depth=F, nperm=1000, carbon.data=Carbondata )#, carbon.age=CarbonFile$age, )

rm(carbon.age)
carbon.age=NULL
by.depth=T
if(by.depth==F | missing(carbon.age)){print(TRUE)} 

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

##########################################################################################################
floco.site.accumulation<-SiteAccumulation("SiteAccumulation")
