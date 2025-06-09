
#load required packages
library(rrBLUP)
library(BMTME)
library(tidyverse)
library(readxl)
library(dplyr)
library(writexl)

#prepare data
data <- as.data.frame(read.csv("Training.csv"))
geno <- as.matrix(data[,-c(1:2)]) #genotype data only
responseVar <- as.vector(data$responseVar)
ids <- as.list(data$ids) #pull genotype IDs

#Calculate and pull marker effects
EBV <-mixed.solve(y=responseVar, Z=geno, K=NULL, SE=FALSE, return.Hinv=FALSE) # 
markerEffects <- EBV$u #u isolates the random effects (marker effects)

#load candidates and update genotypes so test SNPs are identical to training SNPs
candidates <- read.csv("parentCandidates.csv")
candidates = data.frame(candidates[,c(2:3)]) #market class and id only
idCand <- ids[!is.na(ids)]
rownames(candidates) = idCand

rownames(geno) = ids
candidates = merge(candidates,geno, by ="row.names") #match candidates with their genotype

BV <- as.data.frame(data$responseVar) #pull response variable 
BV <- as.data.frame(BV[-nrow(BV),]) #remove unnamed entry at end of DF
rownames(BV) = idCand

#testing only crans
marketclass <- candidates[(candidates$MarketClass == 'cran'),] #pull only crans 
candGeno <- as.matrix(marketclass[,-c(1:3)]) #pull geno only
sampleIDs <- as.data.frame(marketclass$IDS) #pull IDs
GEBVs <- as.data.frame(candGeno %*% markerEffects) #compute GEBVs

#create outupt
GEBVs <- cbind(sampleIDs, GEBVs) #associate GEBVs with IDs
colnames(GEBVs) <- c("sampleIDs", "GEBVs")
rownames(GEBVs) <- GEBVs$sampleIDs
Top <- GEBVs %>% arrange(desc((GEBVs)))
Top[1:5,]

write.csv(Top, "cranGEBVsJan2023.csv")

#calculate accuracy
eval = merge(GEBVs,BV,by="row.names") #merge true and estimated BVs for eval
colnames(eval) = c("id","id","gebv","ebv")
acc = cor(eval$gebv,eval$ebv)
acc

