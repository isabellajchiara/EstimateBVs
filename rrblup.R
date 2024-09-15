
#load required packages
library(rrBLUP)
library(BMTME)
library(tidyverse)
library(readxl)
library(dplyr)
library(writexl)

#prepare data
data <- as.data.frame(read_xlsx("finalTrainingSet.xlsx"))
geno <- data[,-c(1:3)] #genotype data only
responseVar <- as.vector(data$Yield)

ids <- as.list(data$ID) #pull genotype IDs

#pull SNPs with missing calls
geno[geno=="--"] <- NA #no calls become NA
genoVar <- geno %>%  select(where(~ n_distinct(.) > 1)) #remove SNPs will no variance  
M <- as.matrix(genoVar %>%select_if(~!any(is.na(.)))) #remove SNPs with no calls 

#Calculate and pull marker effects
EBV <-mixed.solve(y=responseVar, Z=M, K=NULL, SE=FALSE, return.Hinv=FALSE) # 
markerEffects <- EBV$u #u isolates the random effects (marker effects)

#load candidates and update genotypes so test SNPs are identical to training SNPs
candidates <- read_xlsx("Candidates.xlsx", sheet=2)
candidates = data.frame(candidates[,c(1:2)]) #market class and id only
candidates = candidates[!duplicated(candidates$IDS),] #get individual observations 
candidates = candidates[-nrow(candidates),]
idCand <- ids[!is.na(ids)]
rownames(candidates) = idCand

rownames(M) = ids
M = data.frame(M)
candidates = merge(candidates,M, by ="row.names") #match candidates with their genotype

pheno <- as.data.frame(data$Yield) #pull response variable 
pheno <- pheno[-nrow(pheno),]
pheno <- as.vector(pheno)
BV = as.data.frame(pheno)
rownames(BV) = idCand

#testing only crans
marketclass <- candidates[(candidates$MarketClass == 'cran'),] #pull only crans 
candGeno <- as.matrix(marketclass[,-c(1:3)]) #pull geno only

sampleIDs <- as.data.frame(marketclass$IDS) #pull IDs
GEBVs <- as.data.frame(candGeno %*% markerEffects) #compute GEBVs

GEBVs <- cbind(sampleIDs, GEBVs) #associate GEBVs with IDs
colnames(GEBVs) <- c("sampleIDs", "GEBVs")
rownames(GEBVs) <- GEBVs$sampleIDs
Top <- GEBVs %>% arrange(desc((GEBVs)))
write.csv(Top, "cranGEBVsJan2023.csv")


eval = merge(GEBVs,BV,by="row.names") #merge true and estimated BVs for eval
colnames(eval) = c("id","id","gebv","ebv")
acc = cor(eval$gebv,eval$ebv)
acc

