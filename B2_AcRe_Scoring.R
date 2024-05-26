library(readr)

B2AcReCorellations <- as.data.frame(
  read_csv("~/Desktop/NLME/Null_B2AcRe/B2AcReCorellations.csv"))

# SummaryB2AcRe_Summary<-as.data.frame(
#   read_csv("~/Desktop/NLME/Null_B2AcRe/SummaryB2AcRe_Summary.csv"))
# 
# bestruns<-SummaryB2AcRe_Summary$runs

bestruns<-c("23CC1", "13CC1", "23DC1", "225CC1", "213CC1", "325AC1", " 322AC1", "17CC1", "27CC1", "013EC1", 
            "03DC1", "114BC1", "230CC1", "129BC1", "230DC1")

Scores<-as.data.frame(matrix(data=NA, nrow=length(bestruns), ncol=2))
colnames(Scores)<-c("Run", "Score")

Scores$Run=bestruns

for (i in 1:length(bestruns)){
  currentrun<-B2AcReCorellations[which(B2AcReCorellations$run==bestruns[i]), ]
  
  x<-which(!is.nan(currentrun$rses) & abs(as.numeric(currentrun$rses))>=50)
  
  Scores$Score[i]=sum(currentrun$rses[x])
}

write.csv(Scores, file="ScoresB2AcRe.csv")
