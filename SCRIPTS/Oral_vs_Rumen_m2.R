setwd("/home/alakamp/Documents/PhD/Oral_Rumen_Microbiome") #Path folder with data
library(BGLR)
library(ISLR)
library(CMplot)
library(performance)
library(dplyr)
set.seed(2752) 

##Animals####
#This only getting phenotypes from one datafile only works IF 
##you have checked the same animal has the same records for a trait in both files
Y = read.csv("Sample50120RUMENcount.csv") %>%
  select(., c(2, 7, 4, 16, 18))
X = read.csv("Sample50120ORALcount.csv")%>%
  select(., c(2,4))
A = merge(X, Y, by = "CowID")
rm(Y,X)

nrow(A) - sum(table(A$CowID) == 1) #4 records that are not unique
A$CowID[duplicated(A$CowID)] #4 instances of cow 2106, 2 oral & 2 rumen samples
A = A %>%
  filter(., !CowID == 2106) %>% #121 samples left
  filter(., !Adj.WeanWt == 0) #112 animals
table(A$Sex)

colnames(A) = gsub(".x", "_oral", colnames(A))
colnames(A) = gsub(".y", "_rumen", colnames(A))

###ORAL DATA####
#Read in counts
M <- read.csv("ASVs50120ORALcount.csv")

#Change ASV name columns to rownames of df
ASV_names = M[,1]
M[,1] = NULL
rownames(M) = ASV_names

#Filters by prevalence and makes M a matrix
prevalence_threshold = ncol(M) * 0.05 #prevalent in 10% of animals (change decimal to change %)
##Currently only contains 225 ASV
M = M %>%
  filter(., rowSums(. != 0) >= prevalence_threshold) %>%
  as.matrix()

#Adds 1 to all counts, finds RA for each sample, then log10-transforms all values
##Relative abundance to put all samples on equal standing (even out depth differences)
##Log transform to approx. normal distribution as that's what the model is expecting
M = M + 1
for(col in 1:ncol(M)){
  total = sum(M[,col])
  M[,col] = M[,col] / total
  M[,col] = log10(M[,col])
}

#Microbial matrix, centered and scaled, and ready for BRR
Mc = scale(t(M), center = T, scale = T)

####Adj. Weaning Weight#####
y = A %>%
  select(., c(SampleID_oral, Adj.WeanWt))
cat("Number of animals in analysis is:", nrow(y), "\n")

##Get same IDs as M and filter M to same animals and order as pheno
y$SampleID_oral = paste0("X", y$SampleID_oral)
Mc_working = subset(Mc, rownames(Mc) %in% y$SampleID_oral)
Mc_working = Mc_working[order(match(y$SampleID_oral, rownames(Mc_working))),]

##Setting model as Bayesian Ridge Regression
ETA <- list(MRK=list(X = Mc_working, model = "BRR"))

#Run model
Adj_wwt_model <- BGLR(y = y$Adj.WeanWt, 
                      ETA = ETA, 
                      nIter = 12000, 
                      burnIn = 2000, 
                      saveAt = "Adj_wwt_",  
                      verbose = F)

#Convergence check
varB_trace <- scan("Adj_wwt_ETA_MRK_varB.dat")
varB_mcmc <- mcmc(varB_trace, start = 2001, thin = 5)
plot(varB_mcmc) #Hairy catepillar
summary(varB_mcmc)
acfplot(varB_mcmc) #Should descend to close to 0
effectiveSize(varB_mcmc) #greater than 200-1000
heidel.diag(varB_mcmc) #Should say passed twice
geweke.diag(varB_mcmc) #Should be close to 0, (-1.96 to 1.96 is okay)

# microbial variance
microbial_var = Adj_wwt_model$ETA[[1]]$varB
# residual variance
res_var = Adj_wwt_model$varE
# indiv microbial effects
Adj_wwt_model$ETA[[1]]$b
# total microbiome variance
#sigma_Adj_wwt <- sum(apply(Mc_working, 2, var)) * Adj_wwt_model$ETA[[1]]$varB
sigma_Adj_wwt <- ncol(Mc_working) * microbial_var
# microbiability
lower_var_est = (microbial_var - 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)
upper_var_est = (microbial_var + 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)

cat("m^2 estimate is:\n",
    "Lower:", round(( lower_var_est / (lower_var_est + res_var)),2), "\n",
    "Point:", round(sigma_Adj_wwt  / (sigma_Adj_wwt + res_var), 2), "\n",
    "Upper:", round(( upper_var_est / (upper_var_est + res_var)),2), "\n"
)

####Adj. IMF#####
y = A %>%
  select(., c(SampleID_oral, US.Adj.IMF))
cat("Number of animals in analysis is:", nrow(y), "\n")

##Get same IDs as M and filter M to same animals and order as pheno
y$SampleID_oral = paste0("X", y$SampleID_oral)
Mc_working = subset(Mc, rownames(Mc) %in% y$SampleID_oral)
Mc_working = Mc_working[order(match(y$SampleID_oral, rownames(Mc_working))),]

##Setting model as Bayesian Ridge Regression
ETA <- list(MRK=list(X = Mc_working, model = "BRR"))

#Run model
Adj_wwt_model <- BGLR(y = y$US.Adj.IMF, ETA = ETA, nIter = 12000, burnIn = 2000, verbose = F)

# microbial variance
microbial_var = Adj_wwt_model$ETA[[1]]$varB
# residual variance
res_var = Adj_wwt_model$varE
# indiv microbial effects
Adj_wwt_model$ETA[[1]]$b
# total microbiome variance
#sigma_Adj_wwt <- sum(apply(Mc_working, 2, var)) * Adj_wwt_model$ETA[[1]]$varB
sigma_Adj_wwt <- ncol(Mc_working) * microbial_var
# microbiability
lower_var_est = (microbial_var - 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)
upper_var_est = (microbial_var + 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)

cat("m^2 estimate is:\n",
    "Lower:", round(( lower_var_est / (lower_var_est + res_var)),2), "\n",
    "Point:", round(sigma_Adj_wwt  / (sigma_Adj_wwt + res_var), 2), "\n",
    "Upper:", round(( upper_var_est / (upper_var_est + res_var)),2), "\n"
)

###RUMEN DATA####
#Read in counts
M <- read.csv("ASVs50120RUMENcount.csv")

#Change ASV name columns to rownames of df
ASV_names = M[,1]
M[,1] = NULL
rownames(M) = ASV_names

#Filters by prevalence and makes M a matrix
prevalence_threshold = ncol(M) * 0.05 #prevalent in 5% of animals (change decimal to change %)
##Currently is 1353
M = M %>%
  filter(., rowSums(. != 0) >= prevalence_threshold) %>%
  as.matrix()

#Adds 1 to all counts, finds RA for each sample, then log10-transforms all values
M = M + 1
for(col in 1:ncol(M)){
  total = sum(M[,col])
  M[,col] = M[,col] / total
  M[,col] = log10(M[,col])
}

#Microbial matrix, centered and scaled and ready for BRR
Mc = scale(t(M), center = T, scale = T)

####Adj. Weaning Weight#####
y = A %>%
  select(., c(SampleID_rumen, Adj.WeanWt))
cat("Number of animals in analysis is:", nrow(y), "\n")

##Get same IDs as M and filter M to same animals and order as pheno
y$SampleID_rumen = paste0("X", y$SampleID_rumen)
Mc_working = subset(Mc, rownames(Mc) %in% y$SampleID_rumen)
Mc_working = Mc_working[order(match(y$SampleID_rumen, rownames(Mc_working))),]

##Setting model as Bayesian Ridge Regression
ETA <- list(MRK=list(X = Mc_working, model = "BRR"))

#Run model
Adj_wwt_model <- BGLR(y = y$Adj.WeanWt, ETA = ETA, nIter = 12000, burnIn = 2000, verbose = F)

#Convergence check
varB_trace <- scan("ETA_MRK_varB.dat")
varB_mcmc <- mcmc(varB_trace, start = 2001, thin = 5)
plot(varB_mcmc) #Hairy catepillar
summary(varB_mcmc)
acfplot(varB_mcmc) #Should descend to close to 0
effectiveSize(varB_mcmc) #greater than 200-1000
heidel.diag(varB_mcmc) #Should say passed twice
geweke.diag(varB_mcmc) #Should be close to 0, (-1.96 to 1.96 is okay)

# microbial variance
microbial_var = Adj_wwt_model$ETA[[1]]$varB
# residual variance
res_var = Adj_wwt_model$varE
# indiv microbial effects
Adj_wwt_model$ETA[[1]]$b
# total microbiome variance
#sigma_Adj_wwt <- sum(apply(Mc_working, 2, var)) * Adj_wwt_model$ETA[[1]]$varB
sigma_Adj_wwt <- ncol(Mc_working) * microbial_var
# microbiability
lower_var_est = (microbial_var - 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)
upper_var_est = (microbial_var + 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)

cat("m^2 estimate is:\n",
    "Lower:", round(( lower_var_est / (lower_var_est + res_var)),2), "\n",
    "Point:", round(sigma_Adj_wwt  / (sigma_Adj_wwt + res_var), 2), "\n",
    "Upper:", round(( upper_var_est / (upper_var_est + res_var)),2), "\n"
)

####Adj. IMF#####
y = A %>%
  select(., c(SampleID_rumen, US.Adj.IMF))
cat("Number of animals in analysis is:", nrow(y), "\n")

##Get same IDs as M and filter M to same animals and order as pheno
y$SampleID_rumen = paste0("X", y$SampleID_rumen)
Mc_working = subset(Mc, rownames(Mc) %in% y$SampleID_rumen)
Mc_working = Mc_working[order(match(y$SampleID_rumen, rownames(Mc_working))),]

##Setting model as Bayesian Ridge Regression
ETA <- list(MRK=list(X = Mc_working, model = "BRR"))

#Run model
Adj_wwt_model <- BGLR(y = y$US.Adj.IMF, ETA = ETA, nIter = 12000, burnIn = 2000, verbose = F)

# microbial variance
microbial_var = Adj_wwt_model$ETA[[1]]$varB
# residual variance
res_var = Adj_wwt_model$varE
# indiv microbial effects
Adj_wwt_model$ETA[[1]]$b
# total microbiome variance
#sigma_Adj_wwt <- sum(apply(Mc_working, 2, var)) * Adj_wwt_model$ETA[[1]]$varB
sigma_Adj_wwt <- ncol(Mc_working) * microbial_var
# microbiability
lower_var_est = (microbial_var - 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)
upper_var_est = (microbial_var + 2*Adj_wwt_model$ETA$MRK$SD.varB) * ncol(Mc_working)

cat("m^2 estimate is:\n",
    "Lower:", round(( lower_var_est / (lower_var_est + res_var)),2), "\n",
    "Point:", round(sigma_Adj_wwt  / (sigma_Adj_wwt + res_var), 2), "\n",
    "Upper:", round(( upper_var_est / (upper_var_est + res_var)),2), "\n"
)