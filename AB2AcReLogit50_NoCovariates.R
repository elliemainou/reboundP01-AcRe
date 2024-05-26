# Use Monolix script to fit all possible models of the late treatment group with the timing
# of dynamics (acute vs rebound) 
# Things that I try: fit c, fix c=3,6, 12 +re for acute infection, use dynamics as categorical covariate
# for all parameters
# Also do this by trying different initial conditions???
# How naming of the runs works: Start with "B2AcRe_". Then:
# First digit is for c: 0- fit c, 1- c=3+re, 2- c=6+re, 3=12+re, 4=24+re
# Second digit is for the categorical covariate
# eg B2AcRe_01: fit c, use dynamics as cat cov for delta
# eg B2AcRe_11: fix c=3+re, use dynamics as cat cov for delta

set.seed(12345)
library(ggplot2)

#Install the package and initialize it
##/Applications/MonolixSuite2023R1.app/Contents/Resources/monolixSuite/connectors/lixoftConnectors.tar.gz
#install.packages("/Applications/MonolixSuite2023R1.app/Contents/Resources/monolixSuite/connectors/lixoftConnectors.tar.gz",
#                 repos = NULL, type="source", INSTALL_opts ="--no-multiarch")
library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix", force = TRUE)

#Load the project
loadProject(projectFile="/Users/ejm6007/Desktop/NLME/B2AcReLogit50/B2AcReLogit50_.mlxtran")

scenario <- getScenario()
scenario$tasks = c(populationParameterEstimation = T,
                   conditionalModeEstimation = T,
                   conditionalDistributionSampling = T,
                   standardErrorEstimation=T,
                   logLikelihoodEstimation=T)
scenario$linearization = FALSE
setScenario(scenario)

setwd("/Users/ejm6007/Desktop/NLME/B2AcReLogit50")

tabestimates <- NULL; tabiters <- NULL
init<-NULL ; initprmvalues<-NULL
indiv<-NULL ; indivprmvalues<-NULL
names<-c("A", "B", "C", "D", "E")

# TablePath <- "/storage/work/ejm6007/B2AcReLogit50"
# output_dir_VPC <-"/storage/work/ejm6007/B2AcReLogit50/B2AcReLogit50_VPC"
# output_dir_Fits <- "/storage/work/ejm6007/B2AcReLogit50/B2AcReLogit50_Fits"
# output_dir_PrmDist <-"/storage/work/ejm6007/B2AcReLogit50/B2AcReLogit50_PrmDist"


# ----------------------------------------------------------------------------
# Fit all parameters, do not fix c to different values:
# store the results in tabestimates
# no categorical covariate
# ----------------------------------------------------------------------------
for (j in 1:5){
  loadProject(projectFile="/Users/ejm6007/Desktop/NLME/B2AcReLogit50/B2AcReLogit50_.mlxtran")
  scenario <- getScenario()
  scenario$tasks = c(populationParameterEstimation = T,
                     conditionalModeEstimation = T,
                     conditionalDistributionSampling = T,
                     standardErrorEstimation=T,
                     logLikelihoodEstimation=T)
  scenario$linearization = FALSE
  setScenario(scenario)
  
  #Try initial parameter values that I specified
  if (j==1){
    popParams<-getPopulationParameterInformation()
  }else{
    popParams<-getPopulationParameterInformation()
    #sample new initial estimates uniformly around original values
    popIni<-sapply(popParams$initialValue, function(initialValue){runif(n=1, min=initialValue/2, max=initialValue*2)})
    #log10beta is negative, so I need to switch the bounds
    popIni[4]<-runif(n=1, max=popParams$initialValue[4]/2, min=popParams$initialValue[4]*2)
    #for omega and alpha keep it to 1
    popIni[6:11]<-1
    #set sampled values as new initial estimates
    newPopParams<- popParams
    newPopParams$initialValue<-popIni
    setPopulationParameterInformation(newPopParams)
  }
  
  RunName<-paste(0, names[j], sep = "")
  
  #Save initial parameter values
  init<-merge(rep(RunName, nrow(getPopulationParameterInformation())), getPopulationParameterInformation(), by="row.names")
  init<-init[-1]
  initprmvalues <- rbind(initprmvalues, init)
  
  #run the estimation
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
  
  #Correlations: ANOVA test for signifant correlations
  # yellow: p-value between 0.05-0.1
  # orange: p-value between 0.01-0.05
  # red: p-value less than 0.001
  
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
  estimates$run <- paste("0", names[j], sep = "")
  estimates$Description<-paste("Fit all prms")
  names(estimates)[names(estimates) == "Row.names"] <- "param"
  tabestimates <- rbind(tabestimates, estimates)
  
  #Plots
  #Plot VPC
  m1<-plotVpc(settings = list(xlab="Time (Days)", ylab="Log10 Viral Load"))
  image_name <- paste("VPC_B2AcReLogit50_0", names[j], sep = "")
  # ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
  #        plot = m1, device = "jpeg")
  ggsave(filename = paste0(image_name, ".jpg"), plot = m1, device = "jpeg")
  
  
  #Plot fitted curves
  m2<-plotIndividualFits(settings = list(xlab="Time (Days)", ylab="Log10 Viral Load",
                                         ylim=c(0, 6.5)))
  image_name <- paste("Fits_B2AcReLogit50_0", names[j], sep = "")
  # ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
  #        plot = m2, device = "jpeg")
  ggsave(filename = paste0(image_name, ".jpg"), plot = m2, device = "jpeg")
  
  
  # Plot Distribution of parameter values
  m3 <- plotParametersDistribution()
  image_name <- paste("PrmDist_B2AcReLogit50_0", names[j], sep = "")
  # ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
  #        plot = m3, device = "jpeg")
  ggsave(filename = paste0(image_name, ".jpg"), plot = m3, device = "jpeg")
  
  
  # #save the project
  # file_name<-paste("B2AcReLogit50_0", names[j], ".mlxtran", sep = "")
  # saveProject(projectFile = file_name)
  #save the table
  write.csv(tabestimates, file="B2AcReLogit50NoCovariates.csv")
  write.csv(initprmvalues, file="B2AcReLogit50NoCovariatesInitPrm.csv")
  write.csv(indivprmvalues, file="B2AcReLogit50NoCovariatesIndivPrm.csv")
}



# ----------------------------------------------------------------------------
# fixing c to different values: fix c to 3, 6, 12 per day and run those models
# store the results in tabestimates
# no categorical covariate
# ----------------------------------------------------------------------------
c_values=c(3, 6, 12)

#Try the initial values I have specified
for(i in 1:length(c_values)){
  
  for (j in 1:5){
    loadProject(projectFile="/Users/ejm6007/Desktop/NLME/B2AcReLogit50/B2AcReLogit50_.mlxtran")
    scenario <- getScenario()
    scenario$tasks = c(populationParameterEstimation = T,
                       conditionalModeEstimation = T,
                       conditionalDistributionSampling = T,
                       standardErrorEstimation=T,
                       logLikelihoodEstimation=T)
    scenario$linearization = FALSE
    setScenario(scenario)
    #Try initial parameter values that I specified
    if (j==1){
      popParams<-getPopulationParameterInformation()
      # new value of c
      popIni <- popParams
      popIni$method[which(popParams$name=="c_pop")]="FIXED"
      popIni$initialValue[which(popIni$name=="c_pop")]=as.numeric(c_values[i])
      print(popIni)
      # set values as new initial estimates
      newpopparams <- popIni
      setPopulationParameterInformation(newpopparams)
    }else{
      popParams<-getPopulationParameterInformation()
      #sample new initial estimates uniformly around original values
      popIni<-sapply(popParams$initialValue, function(initialValue){runif(n=1, min=initialValue/2, max=initialValue*2)})
      #log10beta is negative, so I need to switch the bounds
      popIni[4]<-runif(n=1, max=popParams$initialValue[4]/2, min=popParams$initialValue[4]*2)
      #for omega and alpha keep it to 1
      popIni[6:11]<-1
      #set sampled values as new initial estimates
      newPopParams<- popParams
      newPopParams$initialValue<-popIni
      setPopulationParameterInformation(newPopParams)
      
      # new value of c
      popparams<-getPopulationParameterInformation()
      popini <- popparams
      popini$method[which(popParams$name=="c_pop")]="FIXED"
      popini$initialValue[which(popparams$name=="c_pop")]=as.numeric(c_values[i])
      print(popini)
      # set values as new initial estimates
      newpopparams <- popini
      setPopulationParameterInformation(newpopparams)
    }
    RunName<-paste(i, names[j], sep = "")
    
    #Save initial parameter values
    init<-merge(rep(RunName, nrow(getPopulationParameterInformation())), getPopulationParameterInformation(), by="row.names")
    init<-init[-1]
    initprmvalues <- rbind(initprmvalues, init)
    
    #run the estimation
    print(paste("Run", i, names[j]) )
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
    
    #Correlations: ANOVA test for signifant correlations
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
    estimates$run <- paste(i, names[j], sep = "")
    estimates$Description<-paste("c=", c_values[i], "+re")
    names(estimates)[names(estimates) == "Row.names"] <- "param"
    tabestimates <- rbind(tabestimates, estimates)
    
    #Plots
    #Plot VPC
    m1<-plotVpc(settings = list(xlab="Time (Days)", ylab="Log10 Viral Load"))
    image_name <- paste("VPC_B2AcReLogit50_0", names[j], sep = "")
    # ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
    #        plot = m1, device = "jpeg")
    ggsave(filename = paste0(image_name, ".jpg"), plot = m1, device = "jpeg")
    
    #Plot fitted curves
    m2<-plotIndividualFits(settings = list(xlab="Time (Days)", ylab="Log10 Viral Load",
                                           ylim=c(0, 6.5)))
    image_name <- paste("Fits_B2AcReLogit50_0", names[j], sep = "")
    # ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
    #        plot = m2, device = "jpeg")
    ggsave(filename = paste0(image_name, ".jpg"), plot = m2, device = "jpeg")
    
    
    # Plot Distribution of parameter values
    m3 <- plotParametersDistribution()
    image_name <- paste("PrmDist_B2AcReLogit50_0",names[j], sep = "")
    # ggsave(filename = file.path(output_dir_PrmDist, paste0(image_name, ".jpg")),
    #        plot = m3, device = "jpeg")
    ggsave(filename = paste0(image_name, ".jpg"), plot = m3, device = "jpeg")
    
    
    #save the project
    file_name <- paste("B2AcReLogit50_0", names[j],".mlxtran", sep = "")
    saveProject(projectFile = file_name)
    
    #save the table
    write.csv(tabestimates, file="B2AcReLogit50NoCovariates.csv")
    write.csv(initprmvalues, file="B2AcReLogit50NoCovariatesInitPrm.csv")
    write.csv(indivprmvalues, file="B2AcReLogit50NoCovariatesIndivPrm.csv")
  } #loop for j-- init prms 
} #loop for i-- c values