setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")
library(readr)
library(ggplot2)
library(stringr)
library(stringi)

B2AcReLogit50All <- read_csv("B2AcReLogit50All.csv")
Summary_B2AcReLogit50 <- read_csv("Summary_B2AcReLogit50.csv")


#Histogram of AIC of all models
ID<-unique(B2AcReLogit50All$Run)


####Remove the ones where I fit c
x1<-c()

for (i in 1:length(ID)){
  modelname<-ID[i]
  first_char <- substr(modelname, 1, 1)
  if(grepl("0", first_char)){
    x1[i]=i
  }
}
x1<-x1[!is.na(x1)]
ID<-ID[-x1]


####
x<-c()

for(i in 1:length(ID)){
  x[i]<-which(B2AcReLogit50All$Run==ID[i])[1]
}




AIC<-B2AcReLogit50All$AIC[x]

DeltaAIC<-AIC-min(AIC)



hist(DeltaAIC, col="lightblue", xlab=expression(paste(Delta, "AIC")), main="", cex.lab=1.75, cex.axis=1.25, ylim=c(0,150))


#Check values of c for the best model
BestModelsID<-Summary_B2C5C4$Run


B2C5C4CovariatesAbsLRIndivPrmBestModels<-NULL

for (i in 1:length(BestModelsID)){
  current<-B2C5C4CovariatesAbsLRIndivPrm[which(B2C5C4CovariatesAbsLRIndivPrm$x==BestModelsID[i]), ]
  B2C5C4CovariatesAbsLRIndivPrmBestModels<-rbind(B2C5C4CovariatesAbsLRIndivPrmBestModels, current)
}

write.csv(B2C5C4CovariatesAbsLRIndivPrmBestModels, file="B2C5C4CovariatesAbsLRIndivPrmBestModels.csv")
write.csv(B2AcReLogit50All, file="B2AcReLogit50All.csv")

####Create correspondance matrix between name of the run and covariates
names<-c("A", "B", "C", "D", "E")
parameters <- c("delta", "log10pT0", "c", "log10beta", "tstart")
n=length(parameters)
i=1:n
total_sets <- sum(choose(n, i)) 

PrmCombinations<-as.data.frame(matrix(data=NA, nrow=total_sets, ncol=n)) #store prm combinations
PrmCombinations<- as.data.frame(matrix(data=NA, nrow=n, ncol=n))
PrmCombinations[1:n, 1]<-parameters

for (i in 2:n){
  combinations<-t(combn(parameters, i))
  j=1:i
  PrmCombinations[sum(choose(n, j-1)):sum(choose(n, j)), 1:i]=combinations
}

c_values=c(3, 6, 12)

Dynamicsmatrix<-as.data.frame(matrix(data=NA, nrow=nrow(B2AcReLogit50All), ncol=length(parameters)))
cvector<-as.data.frame(matrix(data=NA, nrow=nrow(B2AcReLogit50All), ncol=1))
corr<-as.data.frame(matrix(data=NA, nrow=nrow(B2AcReLogit50All), ncol=2))
colnames(Dynamicsmatrix)=parameters

ID<-unique(B2AcReLogit50All$Run)

for (i in 1:length(ID)){
  currentrows<-which(B2AcReLogit50All$Run==ID[i])
  modelname<-ID[i]
  strlength1<- str_length(modelname)
  sentence<-B2AcReLogit50All$Description[currentrows[1]]
  
  #c value
  x=as.numeric(str_sub(modelname,1,1))
   if(x==0){
    cvector[currentrows, 1]="fit"
  }else if(x==1){
    cvector[currentrows, 1]="3"
  }else if(x==2){
    cvector[currentrows, 1]="6"
  }else if(x==3){
    cvector[currentrows, 1]="12"
  }
  
  #Dynamics  
  if(strlength1==3){
    xDyn=as.numeric(str_sub(modelname,2,2))
  }else if (strlength1==4){
    third_char <- substr(modelname, 3, 3)
    if (grepl("[0-9]", third_char)){
      xDyn=as.numeric(str_sub(modelname,2,3))
    }else{
      xDyn=as.numeric(str_sub(modelname,2,2))
      corr[currentrows, 1]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][1]
      corr[currentrows, 2]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][2]
      }
  }else if (strlength1==5){
    xDyn=as.numeric(str_sub(modelname,2,2))
    corr[currentrows, 1]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][1]
    corr[currentrows, 2]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][2]
  }else if (strlength1==6){
    xDyn=as.numeric(str_sub(modelname,2,3))
    corr[currentrows, 1]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][1]
    corr[currentrows, 2]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][2]
  }else{xDyn=NA}
  
  if(!is.na(xDyn)){
   if("delta" %in% PrmCombinations[xDyn, ]=="TRUE"){
     Dynamicsmatrix[currentrows, 1]<-"Yes"
  }
  if("log10pT0" %in% PrmCombinations[xDyn, ]=="TRUE"){
    Dynamicsmatrix[currentrows, 2]<-"Yes"
  }
  if("c" %in% PrmCombinations[xDyn, ]=="TRUE"){
    Dynamicsmatrix[currentrows, 3]<-"Yes"
  }
  if("log10beta" %in% PrmCombinations[xDyn, ]=="TRUE"){
    Dynamicsmatrix[currentrows, 4]<-"Yes"
  }
  if("tstart" %in% PrmCombinations[xDyn, ]=="TRUE"){
    Dynamicsmatrix[currentrows, 5]<-"Yes"
  }
  }
 
}

B2AcReLogit50AllTotal<-cbind(B2AcReLogit50All, Dynamicsmatrix, cvector, corr)

write.csv(B2AcReLogit50AllTotal, file="B2AcReLogit50AllTotal.csv")

####Create correspondance matrix between name of the run and covariates for the summary matrix
names<-c("A", "B", "C", "D", "E")
parameters <- c("delta", "log10pT0", "c", "log10beta", "tstart")
n=length(parameters)
i=1:n
total_sets <- sum(choose(n, i)) 

PrmCombinations<-as.data.frame(matrix(data=NA, nrow=total_sets, ncol=n)) #store prm combinations
PrmCombinations<- as.data.frame(matrix(data=NA, nrow=n, ncol=n))
PrmCombinations[1:n, 1]<-parameters

for (i in 2:n){
  combinations<-t(combn(parameters, i))
  j=1:i
  PrmCombinations[sum(choose(n, j-1)):sum(choose(n, j)), 1:i]=combinations
}

c_values=c(3, 6, 12)

Dynamicsmatrix<-as.data.frame(matrix(data=NA, nrow=nrow(Summary_B2AcReLogit50), ncol=length(parameters)))
cvector<-as.data.frame(matrix(data=NA, nrow=nrow(Summary_B2AcReLogit50), ncol=1))
corr<-as.data.frame(matrix(data=NA, nrow=nrow(Summary_B2AcReLogit50), ncol=2))
colnames(Dynamicsmatrix)=parameters

ID<-unique(Summary_B2AcReLogit50$Run)

for (i in 1:length(ID)){
  currentrows<-which(Summary_B2AcReLogit50$Run==ID[i])
  modelname<-ID[i]
  strlength1<- str_length(modelname)
  sentence<-Summary_B2AcReLogit50$Description[currentrows[1]]
  
  #c value
  x=as.numeric(str_sub(modelname,1,1))
  if(x==0){
    cvector[currentrows, 1]="fit"
  }else if(x==1){
    cvector[currentrows, 1]="3"
  }else if(x==2){
    cvector[currentrows, 1]="6"
  }else if(x==3){
    cvector[currentrows, 1]="12"
  }
  
  #Dynamics  
  if(strlength1==3){
    xDyn=as.numeric(str_sub(modelname,2,2))
  }else if (strlength1==4){
    third_char <- substr(modelname, 3, 3)
    if (grepl("[0-9]", third_char)){
      xDyn=as.numeric(str_sub(modelname,2,3))
    }else{
      xDyn=as.numeric(str_sub(modelname,2,2))
      corr[currentrows, 1]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][1]
      corr[currentrows, 2]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][2]
    }
  }else if (strlength1==5){
    xDyn=as.numeric(str_sub(modelname,2,2))
    corr[currentrows, 1]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][1]
    corr[currentrows, 2]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][2]
  }else if (strlength1==6){
    xDyn=as.numeric(str_sub(modelname,2,3))
    corr[currentrows, 1]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][1]
    corr[currentrows, 2]=str_extract_all(sentence, "\\b\\w+\\b")[[1]][tail(seq_along(str_extract_all(sentence, "\\b\\w+\\b")[[1]]), 2)][2]
  }else{xDyn=NA}
  
  if(!is.na(xDyn)){
    if("delta" %in% PrmCombinations[xDyn, ]=="TRUE"){
      Dynamicsmatrix[currentrows, 1]<-"Yes"
    }
    if("log10pT0" %in% PrmCombinations[xDyn, ]=="TRUE"){
      Dynamicsmatrix[currentrows, 2]<-"Yes"
    }
    if("c" %in% PrmCombinations[xDyn, ]=="TRUE"){
      Dynamicsmatrix[currentrows, 3]<-"Yes"
    }
    if("log10beta" %in% PrmCombinations[xDyn, ]=="TRUE"){
      Dynamicsmatrix[currentrows, 4]<-"Yes"
    }
    if("tstart" %in% PrmCombinations[xDyn, ]=="TRUE"){
      Dynamicsmatrix[currentrows, 5]<-"Yes"
    }
  }
  
}

Summary_B2AcReLogit50Total<-cbind(Summary_B2AcReLogit50, Dynamicsmatrix, cvector, corr)

write.csv(Summary_B2AcReLogit50Total, file="Summary_B2AcReLogit50Total.csv")
