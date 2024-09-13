#load required packages
library(rrBLUP)
library(BMTME)
library(tidyverse)
library(readxl)


#prepare data
data <- as.data.frame(read_xlsx("trainingData.xlsx", sheet= 2))
data = data[1:269,] #remove NAs at the bottom of the DF
geno <- data[,-c(1:15)] #genotype data only
pheno <- as.vector(data$Yield) #pull response variable 
ids <- as.list(data$ID) #pull genotype IDs

#pull SNPs with missing calls
geno[geno=="--"] <- NA #no calls become NA
genoVar <- geno %>%  select(where(~ n_distinct(.) > 1)) #remove SNPs will no variance  
M <- as.matrix(genoVar %>%select_if(~!any(is.na(.)))) #remove SNPs with no calls 

#Calculate and pull marker effects
EBV <-mixed.solve(y=pheno, Z=M, K=NULL, SE=FALSE, return.Hinv=FALSE) # 
markerEffects <- EBV$u #u isolates the random effects (marker effects)

#load candidates and update genotypes so test SNPs are identical to training SNPs
candidates <- read_xlsx("Candidates.xlsx", sheet=2)
candidates = data.frame(candidates[,c(1:2)]) #market class and id only
rownames(M) = ids
M = data.frame(M)
candidates = merge(candidates,M, by ="row.names") #match candidates with their genotype

#testing only crans
marketclass <- candidates[(candidates$MarketClass == 'cran'),] #pull only crans 
candGeno <- as.matrix(marketclass[,-c(1:3)]) #pull geno only

sampleIDs <- as.data.frame(marketclass[,3]) #pull IDs
GEBVs <- as.data.frame(candGeno %*% markerEffects) #compute GEBVs

GEBVs <- cbind(sampleIDs, GEBVs) #associate GEBVs with IDs
colnames(GEBVs) <- c("sampleIDs", "GEBVs")
Top <- GEBVs %>% arrange(desc((GEBVs)))
Top

write.csv(Top, "cranGEBVsJan2023.csv")
