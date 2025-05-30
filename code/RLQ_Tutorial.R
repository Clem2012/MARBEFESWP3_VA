#' ---
#' title: "RLQ analysis"
#' author: "Clement Garcia, Cefas"
#' date: "20 April 2025"
#' output:
#'   html_document: default
#'   word_document: default
#'   pdf_document: default
#' ---
#' 
## ----setup, include=FALSE----------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

#' This document provides an example data and R code for performing **RLQ analysis** (RLQ) and **Fourth Corner**, (4thC) RLQ combines three tables: site x species, site x environment and species x trait in one single ordination and 4thC tests the significance of different elements of the three tables.
#' For a complete working example, please refer to the [tutorial from Dray et al. (2016)](https://wiley.figshare.com/articles/dataset/Supplement_1_A_tutorial_to_perform_fourth-corner_and_RLQ_analyses_in_R_/suppl-1.pdf).
#' 
#' # 1. Preparation
#' 
#' ## 1.1. Packages and data
#' 
#' The analyses require the R packages `ade4` and `tidyverse`.
#' 
## ---- message=FALSE----------------------------------------------------------------------------------------------
# Clear memory
rm(list=ls())
# Load packages
library(ade4) ## Perform RLQ & 4thC
library(tidyverse) ## data formatting and plotting (ggplot2)

#' 
#' If you get an error message, check that the R packages are installed correctly. If not, use the command: `install.packages(c("ade4", "tidyverse"))`.
#' 
#' The example dataset is available as the Rdata file `NEAtl_FishTraitEnv.Rdata`, available for download [here](https://github.com/rfrelat/TraitEnvironment/raw/main/NEAtl_FishTraitEnv.Rdata).  
#' 
#' ## 1.2 Load and inspect data
#' 
#' Make sure the file `NEAtl_FishTraitEnv.Rdata` is in your working directory, then load it in R.
#' 
## ----------------------------------------------------------------------------------------------------------------
load("NEAtl_FishTraitEnv.Rdata")

#' 
#' The Rdata file contains four objects: 
#' 
#' - `abu` Species abundance by grid cell (the community matrix)
#' - `env` Environmental conditions by grid cell
#' - `trait` Traits by species
#' - `coo` The spatial coordinates (lat/lon) for each grid cell
#' 
## ----------------------------------------------------------------------------------------------------------------
# Inspect data tables 
dim(abu) 

dim(env)
names(env)

dim(trait)
names(trait)

#' 
#' The `trait` table contains `r ncol(trait)`  traits (i.e variable, in column) characterizing `r nrow(trait)` taxa (in rows). The `r ncol(trait)` traits broadly represent the life history and ecology of fish in terms of their feeding, growth, survival and reproduction. These are:  
#' 
#' - Trophic level
#' - K: the growth rate, calculated as Von Bertalanffy growth coefficient in year$^{-1}$
#' - Lmax: maximum body length in cm
#' - Lifespan
#' - Offspring.size_log: egg diameter, length of egg case or length of pup in mm
#' - Fecundity_log: number of offspring produced by a female per year
#' - Age.maturity: in years
#' 
#' Trait values for fecundity and offspring size were log-transformed to reduce the influence of outliers. 
#' 
#' The `env` table contains `r ncol(env)`  environmental variables (in column) characterizing `r nrow(env)` grid cells (in rows). The environmental variables measure hydrography, habitat, food availability and anthropogenic pressures, which are known to affect the distribution of fish species. These are:  
#' 
#' - Depth: depth in meter, directly measured during the survey.
#' - SBT: monthly sea bottom temperature in °C from the Global Ocean Physics Reanalysis (GLORYSs2v4) 
#' - SBS: monthly sea bottom salinity from the Global Ocean Physics Reanalysis (GLORYSs2v4) 
#' - Chl: Chlorophyll a concentration (in $mg.m^{-3}$) as a proxy for primary production and food availability from the GlobColour database
#' - SBT_sea: seasonality of sea bottom temperature, calculated as the difference between the warmest and the coldest month of the year.
#' - Chl_sea: seasonality of chlorophyll a concentration, calculated as the difference between the highest and the lowest primary production in the year
#' - Fishing: the cumulative demersal fishing pressure in 2013, estimated globally by Halpern et al. 2015, [DOI 10.1038/ncomms8615](https://doi.org/10.1038/ncomms8615). 
#' 
#' Importantly, the rows in `abu` correspond to the same grid cell than the rows in `env`, and the column in `abu` correspond to the same taxa than the rows in `trait`.  If you want to learn how to create such dataset, see the short [tutorial on setting trait-environement dataset](https://rfrelat.github.io/CleanDataTER.html).  
#' 
## ----------------------------------------------------------------------------------------------------------------
all(row.names(abu)==row.names(env))
all(colnames(abu)==row.names(trait))

#' 
#' Using this fish community of the Northeast Atlantic as an example, you will learn how to estimate CWM traits and diversity indices, as well as plotting patterns and identifying key drivers and state-pressure relationships (using GAMs and random forest)
#' 
#' 
#' # 2. Calculate Community-weighted mean (CWM) traits and various diversity indices 
#' 
#' ## 2.1 Calculate diversity indices
#' 
#' We use the function `FD()` from the `FD package`. Please see the help `?dbFD` to inspect the function and all its options. The option `stand.FRic=T` means that FRic will range between 0 and 1
#' 
#' Be aware, the computation of `dbFD` can take several minutes.
## ----------------------------------------------------------------------------------------------------------------
fishFD<-dbFD(trait, abu, stand.FRic=T)  

#' 
#' ## 2.2 Inspect results
#' 
## ----------------------------------------------------------------------------------------------------------------
names(fishFD)

#' 
#' 
#' - nbsp is the Species Richness              
#' - FRic is the Functional Richness           
#' - FEve is the Functional Eveness            
#' - FDiv is the Functional Divergence         
#' - RaoQ is the Functional Entropy            
#' - CWM contains the Community weighted mean traits 
#' 
#' ## 2.3 Reformat the output
#' 
## ----------------------------------------------------------------------------------------------------------------
# Reformat into an output data frame
fishFD<-do.call(cbind, fishFD)

## Add lon, lat to output data frame
fishFD[,c("lon","lat")]<-coo

#' 
#' 
#' ## 2.4 Map biodiversity indices
#' 
## ----------------------------------------------------------------------------------------------------------------
colpal<-rev(brewer.pal(11,"RdYlBu")) # Color palette for the map
#We plot the CWM of Trophic level
mapggplot(fishFD$lon, fishFD$lat, 
          fishFD$CWM.Trophic.level, colpal = colpal)


#' 
#' # 3. Statistically investigate trait-environment relationships
#' 
#' Here, we model with GAMS and random forest the trait responses at the community level (using CWM).
#' 
#' 
#' ## 3.1 Preparation of the dataset
#' 
#' We want to model CWM from environmental data, so we need to merge these information in one data frame.
#' 
## ----------------------------------------------------------------------------------------------------------------
# Merge fishFD with env data frame
# Because the rows of fishFD match the rows of env, we can combine them
FDenv<-cbind(fishFD,env)
head(FDenv)

#' ## 3.2 Simple exploratory GAMs
#' 
#' We use the function `gam()` from the `mgcv package`. Please see the documentation `?gam` to inspect the function and all its options. 
#' 
#' 
#' We define a function to quickly plot GAM partial smooth plots and return summary statistics for the fitted model
#' 
#' 
## ----------------------------------------------------------------------------------------------------------------
funGam<-function(Var){
  colnames(FDenv)[which(colnames(FDenv)==Var)]<-"Var"
  fit<-gam(Var~s(Depth, k=3)+s(SBT,k=3)+s(SBS,k=3)+s(Chl,k=3)+s(SBT_sea,k=3)+ s(Chl_sea, k=3)+s(Fishing, k=3), na.action=na.exclude, data=FDenv)
  par(mfrow=c(4,2),mar=c(4,4,0.2,0.2))
  plot(fit,shade=T,shade.col="grey",res=T,rug=F,pch=20,ylab="")
  par(mfrow=c(1,1),mar=c(5,4,2,2))
  return(summary(fit))
}

#' 
#' Now, let's try to model trophic level from environmental data 
#' 
## ----------------------------------------------------------------------------------------------------------------
# list all output of FD calculation
Metrics<-colnames(fishFD)
Metrics

#Plot GAM partial smooth plots and return summary statistics for the fitted model
funGam(Metrics[9])

#' 
#' ## 3.3 Detailed and thorough model selection
#' 
#' ### 3.3.1 Considering potential spatial autocorrelation
#' 
#' Here we use a mixed GAM approach taking into account potential spatial autocorrelation (by allowing for correlated error structures)
#' We use the function `uGamm()` from the `MuMIn package` which is a wrapper function for gamm. Please see the documentation `?uGamm` to inspect the function and all its options. 
#' 
#' We will use trophic level as an example. First, we check distribution of response (below with `hist()`). Consider log transforming or change family statement in gam in case of strong deviations from normal distribution.
#' 
## ----------------------------------------------------------------------------------------------------------------
# Check distribution of response
hist(FDenv$CWM.Trophic.level) 

#' 
#' Then we can compute the mixed GAM considering potential spatial autocorrelation.  
#' 
#' **Be aware, the computation can take several minutes.**
#' 
## ----eval=FALSE--------------------------------------------------------------------------------------------------
## # Create dummy variable to include as random factor
dummy <- rep(1, dim(FDenv)[1])
## ### Dummy variables to be add as random factor (mandatory, but won't 'change' anything)
## 
fitTL <- uGamm(CWM.Trophic.level ~ s(Depth, k=3)
                    +s(SBT,k=3)
                    +s(SBS,k=3)
                    +s(Chl,k=3)
                    +s(SBT_sea,k=3)
                    +s(Chl_sea, k=3)
                    +s(Fishing, k=3),
                    gaussian(link = "identity"),random = list(dummy=~1), correlation = corGaus(form = ~ lon+lat), # Can try other error structures
                    data = FDenv, control=lmeControl(opt="optim"))
## ----------------------------------------------------------------------------------------------------------------
# Inspect summary statistics
summary(fitTL$lme)
summary(fitTL$gam)


#' 
#' 
## ----eval=TRUE---------------------------------------------------------------------------------------------------
# Check diagnostics
par(mfrow=c(2,2))
gam.check(fitTL$gam)

#' 
## ----eval=TRUE---------------------------------------------------------------------------------------------------
par(mfrow=c(4,2),mar=c(4,4,0.2,0.2))
plot(fitTL$gam,shade=T,shade.col="grey",res=T,rug=F,pch=20)
par(mfrow=c(1,1),mar=c(5,4,2,2))

#' 
#' ### 3.3.2 Automated model selection
#' 
#' Here, we perform automated model selection testing all combinations of predictors using the `dredge()` function (and rank according to AICc)
#' Please see the documentation `?dredge` from the `MuMIn package` to inspect the function and all its options.
#' 
#' **Be aware, the computation can take quite some time.**
#' 
## ----eval=FALSE--------------------------------------------------------------------------------------------------
#results <- dredge(fitTL, m.lim=c(1,4), rank="AICc") # Here test max 4 variables per model to reduce run time

## ----echo=FALSE--------------------------------------------------------------------------------------------------
# Or just load and view the results
load("results.Rdata")
results

## ----------------------------------------------------------------------------------------------------------------
subset(results, delta <5)  # Depth, Fishing, SBT, SBT_sea seem to be key variables
# Calculate and view relative variable importance (RVI) scores
importance(results) 

#' 
## ----eval=FALSE--------------------------------------------------------------------------------------------------
# Fit and inspect the "best" model
fitTLb <- uGamm(CWM.Trophic.level ~ s(Depth, k=3)
                +s(Fishing,k=3)
                +s(SBT,k=3)
                +s(SBT_sea,k=3),
                gaussian(link = "identity"),random = list(dummy=~1), correlation = corGaus(form = ~ lon+lat), # Can try other error structures
                data = FDenv, control=lmeControl(opt="optim"))
## 
## ----------------------------------------------------------------------------------------------------------------
# Inspect summary statistics
summary(fitTLb$lme)
summary(fitTLb$gam)

#' 
## ----eval=TRUE---------------------------------------------------------------------------------------------------
# Check diagnostics
par(mfrow=c(2,2))
gam.check(fitTLb$gam)

#' 
## ----eval=TRUE---------------------------------------------------------------------------------------------------
# Plot model smooth terms
par(mfrow=c(2,2),mar=c(4,4,0.2,0.2))
plot(fitTLb$gam,shade=T,shade.col="grey",res=T,rug=F,pch=20)

#' 
#' 
#' # 4. Random forest
#' 
#' Instead of GAM, we can try to model the CWM trait with random forest. 
#' We use the function `randomForest()` from the `randomForest package`. Please see the documentation `?randomForest` to inspect the function and all its options. 
#' 
#' Again we create a function to compute random forests and plot the outputs.
## ----------------------------------------------------------------------------------------------------------------
# Wrapper function to explore RF for a given trait
funRf<-function(Var){
  colnames(FDenv)[which(colnames(FDenv)==Var)]<-"Var"
  fit<-randomForest(Var~Depth+SBT+SBS+Chl+SBT_sea+Chl_sea+Fishing, data=FDenv,ntree=1000,importance=T, mtry=2)
  par(mfrow=c(4,2),mar=c(4,4,0.2,0.2))
  partialPlot(fit, x.var=Depth,FDenv,main="")
  partialPlot(fit, x.var=SBT,FDenv,main="")
  partialPlot(fit, x.var=SBS,FDenv,main="")
  partialPlot(fit, x.var=Chl,FDenv,main="")
  partialPlot(fit, x.var=SBT_sea,FDenv,main="")
  partialPlot(fit, x.var=Chl_sea,FDenv,main="")
  partialPlot(fit, x.var=Fishing,FDenv,main="")
  par(mfrow=c(1,1),mar=c(5,4,2,2))
  return(fit)
}

#' 
#' Now, we compute the randomForest for CWM trophic level.
## ----------------------------------------------------------------------------------------------------------------
fitTLrf <- funRf(Metrics[9])# Plot random forest response plots and summary stats

#' 
#' #### Check error against number of trees used for training
## ----------------------------------------------------------------------------------------------------------------
plot(fitTLrf) 

#' It is stable after 200 (hence 1000 trees is more than enough).
#' 
#' #### Inspect summary stats
## ----------------------------------------------------------------------------------------------------------------
print(fitTLrf) # Inspect summary stats

#' 
#' #### Check variable importance plots
#' 
## ----------------------------------------------------------------------------------------------------------------
varImpPlot(fitTLrf)

########################################################### 
##########################################################
#' # Optional - Generalized linear mixed models (GLMM)
##########################################################
#' We test GLMM and extract parameters to facilitate comparisons across organism groups and areas. 
#' 
#' Please note that you may need to introduce square terms to account for non-linearities.
#' 
#' Additionally, you need to load the package `nlme`.
#' 
## ----------------------------------------------------------------------------------------------------------------
library(nlme)

#' 
## ---- eval=FALSE-------------------------------------------------------------------------------------------------
 fitTLglm <- lme(CWM.Trophic.level~Depth+SBT+SBS+Chl+SBT_sea+Chl_sea+Fishing,
                 random = list(dummy=~1), correlation = corGaus(form = ~ lon+lat),
                 data = FDenv, control=lmeControl(opt="optim"))

#----------------------------------------------------------------------------------------------------------------
summary(fitTLglm)
anova(fitTLglm)

# Extract parameters
summary(fitTLglm)

summary(fitTLglm)$coefficients

##Plot some diagnostics
plot(fitTLglm)
qqnorm(residuals(fitTLglm))

#observed versus fitted values
plot(fitTLglm, CWM.Trophic.level~ fitted(.), abline = c(0,1))

#########################
# END
#################################################
