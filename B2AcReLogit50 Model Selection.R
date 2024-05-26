setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")
library(readr)
library(readxl)
B2AcReLogit50All <- read_csv("B2AcReLogit50AllTotal.csv")


#How many models I tried 
initial<-unique(B2AcReLogit50All$Run)
length(initial)


# 1. Remove runs where there is an NaN in the RSE-- 
remove1<-unique(B2AcReLogit50All$Run[which(is.na(B2AcReLogit50All$rses))])
length(remove1)

remaining1<-setdiff (initial, remove1)

B2AcReLogit50All1 <- subset(B2AcReLogit50All, Run %in% remaining1)

# 2. Keep those with AIC up to 10 points higher of min AIC

AIC_threshold=min(B2AcReLogit50All1$AIC)+10


B2AcReLogit50All2 <-B2AcReLogit50All1[which(B2AcReLogit50All1$AIC<=AIC_threshold),  ]
length(unique(B2AcReLogit50All2$Run))


write.csv(B2AcReLogit50All2, file="B2AcReLogit50AllTotal2.csv")
ID<-unique(B2AcReLogit50All2$Run)


ID<-ID[order(ID, decreasing=FALSE)]
ID

#Now look at the VPC and remove the ones with problematic VPC

#Create summary table 

B2AcReLogit50AllTotal2 <- read_csv("B2AcReLogit50AllTotal2.csv")


ID<-unique(B2AcReLogit50AllTotal2$Run)
names<-c("Run", "AIC", "Red", "Orange", "Yellow", "Score", "Description", 
         "delta",	"log10pT0",	"c",	"log10beta",	"tstart",	"cvalue",	"corr1",	"corr2")
Summary_B2AcReLogit50<-as.data.frame(matrix(data=NA, nrow=length(ID), ncol=length(names)))
colnames(Summary_B2AcReLogit50)<-names

#Add the scoring system:sum up all the absolute values of the RSEs and divide
#by the number of parameters in the model

for (i in 1:length(ID)){
  currentmodel<-B2AcReLogit50AllTotal2[which(B2AcReLogit50AllTotal2$Run==ID[i]), ]
 
   Summary_B2AcReLogit50$Run[i]<-currentmodel$Run[1]
   Summary_B2AcReLogit50$AIC[i]<-currentmodel$AIC[1]
   Summary_B2AcReLogit50$Red[i]<-currentmodel$rsesColorCounts[1]
   Summary_B2AcReLogit50$Orange[i]<-currentmodel$rsesColorCounts[2]
   Summary_B2AcReLogit50$Yellow[i]<-currentmodel$rsesColorCounts[3]
   Summary_B2AcReLogit50$Description[i]<-currentmodel$Description[1]
   x<-which(!is.nan(currentmodel$rses) & abs(as.numeric(currentmodel$rses))>=50)
   Summary_B2AcReLogit50$Score[i]<-sum(abs(as.numeric(currentmodel$rses)))/nrow(currentmodel)
   Summary_B2AcReLogit50$delta[i]<-currentmodel$delta[1]
   Summary_B2AcReLogit50$log10pT0[i]<-currentmodel$log10pT0[1]
   Summary_B2AcReLogit50$c[i]<-currentmodel$c[1]
   Summary_B2AcReLogit50$log10beta[i]<-currentmodel$log10beta[1]
   Summary_B2AcReLogit50$tstart[i]<-currentmodel$tstart[1]
   Summary_B2AcReLogit50$cvalue[i]<-currentmodel$cvalue[1]
   Summary_B2AcReLogit50$corr1[i]<-currentmodel$corr1[1]
   Summary_B2AcReLogit50$corr2[i]<-currentmodel$corr2[1]
   
}

write.csv(Summary_B2AcReLogit50, file="Summary_B2AcReLogit50New1.csv")





