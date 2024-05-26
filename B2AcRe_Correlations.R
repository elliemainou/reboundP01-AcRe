module use /storage/icds/RISE/sw8/modules/
module load r
R

set.seed(12345)
library(ggplot2)
library(lixoftConnectors)
library(readr)
library(stringr)

#Install the package and initialize it
initializeLixoftConnectors(software = "monolix", force = TRUE)

tabestimates <- NULL; tabiters <- NULL
init<-NULL ; initprmvalues<-NULL
indiv<-NULL ; indivprmvalues<-NULL
names<-c("A", "B", "C", "D", "E")

# We want to test a categorical covariate on all possible combinations of parameters.
# Here we create a matrix containing all these combinations (n choose k)

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

TablePath <- "/storage/work/ejm6007/B2AcRe/B2AcRe_Tables"
output_dir_VPC <-"/storage/work/ejm6007/B2AcRe/B2AcRe_VPC"
output_dir_Fits <- "//storage/work/ejm6007/B2AcRe/B2AcRe_Fits"
output_dir_PrmDist <-"/storage/work/ejm6007/B2AcRe/B2AcRe_PrmDist"

#Load matrices
B2AcRe <- read_csv("/storage/work/ejm6007/B2AcRe/B2AcRe_Tables/B2AcReNoCorr.csv")
B2AcReInitPrm <-  read_csv("/storage/work/ejm6007/B2AcRe/B2AcRe_Tables/B2AcReNoCorrInitPrm.csv")

# B2AcRe <- read_csv("Downloads/B2AcReNoCorr.csv")
# B2AcReInitPrm <-  read_csv("Downloads/B2AcReNoCorrInitPrm.csv")

#Find for which runs I need to add correlations
B2AcReToAddCorr<-B2AcRe[which(B2AcRe$CorrelationsRed!= "NA - NA"),]

ID<-B2AcReToAddCorr$Run

for (z in 101:length(ID)){
  #Select the initial prm values
  initprms<-B2AcReInitPrm[which(B2AcReInitPrm$x==ID[z]), ]
  modelresults<-B2AcRe[which(B2AcRe$Run==ID[z]), ]
  
  corrcurrent<-modelresults$CorrelationsRed[which(!is.na(modelresults$CorrelationsRed))]
  
  for(l in 1:length(corrcurrent)){
    #which will be the categorical covariates
    strlength1<- str_length(modelresults$Run[1])
    
    k= as.numeric(str_sub(ID[z], 2,(strlength1-1)))
    i=as.numeric(str_sub(ID[z], 1,1))
    j=which(names==str_sub(ID[z], strlength1,strlength1))
    
    #load project
    loadProject(projectFile="/storage/work/ejm6007/B2AcRe/B2AcRe_3.mlxtran")
    #loadProject(projectFile="/Users/ejm6007/Desktop/NLME/B2AcRe/B2AcRe_.mlxtran")
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    
    #specify initial parameter values
    popParams<-getPopulationParameterInformation()
    #sample new initial estimates uniformly around original values
    popIni<-c()
    for (m in 1:5){ 
      popIni[m]<-initprms$initialValue[which(initprms$name==popParams$name[m])]
    }#loop for m
    #for omega and alpha keep it to 1
    popIni[6:11]<-1
    #set sampled values as new initial estimates
    newPopParams<- popParams
    newPopParams$initialValue<-popIni
    setPopulationParameterInformation(newPopParams)
    popparams<-getPopulationParameterInformation()
    popini <- popparams
    if ( i>0){# fix value of c
      popini$method[which(popParams$name=="c_pop")]="FIXED"
    }
    print(popini)
    # set values as new initial estimates
    newpopparams <- popini
    setPopulationParameterInformation(newpopparams)
    
    #### Add categorical covariate for Abs   
    if (!is.na(k)){
      #Add categorical covariates
      currentcov<-!is.na(PrmCombinations[k, ])
      
      x<-as.character(PrmCombinations[k, which(!is.na(PrmCombinations[k, ]))]) #parameters that will be tested
      number<-length(which(!is.na(PrmCombinations[k, ])))
      # Create a data frame with parameter names and values
      df <- data.frame(x, rep("Dynamics", number), rep(TRUE, number))
      colnames(df) <- c("names", "Covariate" ,"values")
      
      # Use the data frame to generate a string with the new values
      args <- ""
      for (row in 1:nrow(df)) {
        args <- paste(args, df$names[row], "= c(",  df$Covariate[row], "=", df$values[row], "),")
      }
      # Remove the last comma
      args <- substr(args, 1, nchar(args) - 1)
      # Wrap the args in the function call
      args <- paste("setCovariateModel(", args, ")")
      # Print what we have as a sanity check and attempt to parse it
      #print(args)
      eval(parse(text = args))
    }
    #What is the current iteration 
    RunName<- paste(i, k, names[j], "C", l, sep = "")
    
    #Save initial parameter values
    init<-merge(rep(RunName, nrow(getPopulationParameterInformation())), getPopulationParameterInformation(), by="row.names")
    init<-init[-1]
    initprmvalues <- rbind(initprmvalues, init)
    
    #### Add correlations
    corr1<-word(corrcurrent[1], 1)
    corr2<-word(corrcurrent[1], 3)
    
    setCorrelationBlocks(id = list( c(corr1,corr2) ) )
    
    # run the estimation
    print(RunName)
    runScenario()
    
    #Save individual parameter values
    indiv<-merge(rep(RunName, nrow(getEstimatedIndividualParameters()$conditionalMode)), 
                 getEstimatedIndividualParameters()$conditionalMode, by="row.names")
    indiv<-indiv[-1]
    indivprmvalues <- rbind(indivprmvalues, indiv)
    
    # store the estimates and s.e. in a table
    #parameter estimates
    estimates <- as.data.frame(getEstimatedPopulationParameters())
    names(estimates) <- "estimate"
    
    #rise in standard error
    rses <- getEstimatedStandardErrors()$stochasticApproximation$rse
    names(rses) <- getEstimatedStandardErrors()$stochasticApproximation$parameter
    rses <- as.data.frame(rses)
    estimates <- merge(estimates, rses, by="row.names")
    
    #RSE, which are the unstable values
    estimates$rsesColor<-"NA"
    estimates$rsesColor[which(!is.nan(estimates$rses) & abs(as.numeric(estimates$rses))>=50 
                              & abs(as.numeric(estimates$rses))<100)]="Yellow" 
    estimates$rsesColor[which(!is.nan(estimates$rses) & abs(as.numeric(estimates$rses))>=100 
                              & abs(as.numeric(estimates$rses))<200)]="Orange" 
    estimates$rsesColor[which(!is.nan(estimates$rses) & abs(as.numeric(estimates$rses))>=200)]="Red" 
    
    #count them
    estimates$rsesColorCounts=NA
    estimates$rsesColorCounts[1]=paste("Red=", length(which(estimates$rsesColor=="Red")))
    estimates$rsesColorCounts[2]=paste("Orange=", length(which(estimates$rsesColor=="Orange")))
    estimates$rsesColorCounts[3]=paste("Yellow=", length(which(estimates$rsesColor=="Yellow")))
    
    #AIC
    estimates$AIC <- round(getEstimatedLogLikelihood()$importanceSampling[2], 2)
    
    #Correlations: ANOVA test for significant correlations
    # yellow: p-value between 0.05-0.1
    # orange: p-value between 0.01-0.05
    # red: p-value less than 0.001
    
    #Which sets of parameters have significant correlations?
    
    #Yellow
    x_yellow=which(getTests()$randomEffectsCorrelation$p.value>0.05 & getTests()$randomEffectsCorrelation$p.value<0.1)
    estimates$CorrelationsYellow=NA
    for(y in 1:length(x_yellow)){
      estimates$CorrelationsYellow[y]=paste(getTests()$randomEffectsCorrelation$eta1[x_yellow[y]], "-", getTests()$randomEffectsCorrelation$eta2[x_yellow[y]])
      estimates$CorrelationsYellow[y]=gsub("eta_","",as.character(estimates$CorrelationsYellow[y]))
    }
    
    #Orange
    x_orange=which(getTests()$randomEffectsCorrelation$p.value>0.01 & getTests()$randomEffectsCorrelation$p.value<0.05)
    estimates$CorrelationsOrange=NA
    
    for(o in 1:length(x_orange)){
      estimates$CorrelationsOrange[o]=paste(getTests()$randomEffectsCorrelation$eta1[x_orange[o]], "-", getTests()$randomEffectsCorrelation$eta2[x_orange[o]])
      estimates$CorrelationsOrange[o]=gsub("eta_","",as.character(estimates$CorrelationsOrange[o]))
    }
    
    #Red
    x_red=which(getTests()$randomEffectsCorrelation$p.value<0.01)
    estimates$CorrelationsRed=NA
    for(r in 1:length(x_red)){
      estimates$CorrelationsRed[r]=paste(getTests()$randomEffectsCorrelation$eta1[x_red[r]], "-", getTests()$randomEffectsCorrelation$eta2[x_red[r]])
      estimates$CorrelationsRed[r]=gsub("eta_","",as.character(estimates$CorrelationsRed[r]))
    }
    
    #Description of the model
    estimates$Run <-RunName
    if(i==0){
      text1<-paste("c=fit, ")
    }else {text1<-paste("c=", c_values[i], "+re, ", sep = "")}
    
    if(!is.na(k)){
      text2<-paste("Dynamics= cat cov ", paste(x, collapse = ","))
    }else{text2<-paste("no cov")}
    
    estimates$Description<-paste(text1, text2, ", ",
                                 names[j], ", corr: ", paste(corr1), "-", paste(corr2), sep = "")
    print(paste(text1, text2, ", ",
                names[j], ", corr: ", paste(corr1), "-", paste(corr2), sep = ""))
    names(estimates)[names(estimates) == "Row.names"] <- "param"
    tabestimates <- rbind(tabestimates, estimates)
    
    #Plots
    #Plot VPC
    computeChartsData(exportVPCSimulations = TRUE)
    m1<-plotVpc(settings = list(xlab="Time (Days)", ylab="Log10 Viral Load"))
    image_name <- paste("VPC_B2AcRe_", i, k, names[j],"C", l, sep = "")
    ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
           plot = m1, device = "jpeg")
    
    #Plot fitted curves
    m2<-plotIndividualFits(settings = list(xlab="Time (Days)", ylab="Log10 Viral Load",
                                           ylim=c(0, 6.5)))
    image_name <- paste("Fits_B2AcRe_", i, k, names[j], "C", l, sep = "")
    ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
           plot = m2, device = "jpeg")
    
    # Plot Distribution of parameter values
    m3 <- plotParametersDistribution()
    image_name <- paste("PrmDist_B2AcRe_", i, k, names[j],"C", l, sep = "")
    ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
           plot = m3, device = "jpeg")
    
    #save the table
    write.csv(tabestimates,  file.path(TablePath,"B2AcReCovariatesAbsLRCorr3.csv"))
    write.csv(initprmvalues,  file.path(TablePath,"B2AcReCovariatesAbsLRInitPrmCorr3.csv"))
    write.csv(indivprmvalues,  file.path(TablePath,"B2AcReCovariatesAbsLRIndivPrmCorr3.csv"))
    
  }#loop for l-- number of red correlations
} #loop for z-- IDs  
}# loop for k -- which is the current r