#load required packages
library(rrBLUP)
library(BMTME)
library(tidyverse)
library(readxl)


#prepare data
data <- as.data.frame(read_xlsx("trainingData.xlsx", sheet= 2))
geno <- data[,-c(1,2)]
pheno <- as.data.frame(data$Yield)
pheno <- as.matrix(pheno[1:269,]) #remove NAs at the bottom of the DF
ids <- as.data.frame(data$ID)
ids = ids[1:269,] #remove NAs at the bottom of the DF

#pull SNPs with missing calls
geno[geno=="--"] <- NA
genoVar <- geno %>%  select(where(~ n_distinct(.) > 1))
genoVar <- genoVar[1:269,]
M <- as.matrix(genoVar %>%select_if(~!any(is.na(.))))

#Calculate and pull marker effects
EBV <-mixed.solve(y=pheno, Z=M, K=NULL, SE=FALSE, return.Hinv=FALSE)
markerEffects <- EBV$u #u isolates the coefficients (marker effects)

#load candidates and update genotypes so test SNPs are identical to training SNPs
candidates <- read_xlsx("Candidates.xlsx", sheet=2)
candidates = candidates[,c(1:2)] #market class and id only
candidates = data.frame(candidates)
idCand = candidates$IDS
rownames(candidates) = idCand
rownames(M) = ids
M = data.frame(M)
candidates = merge(candidates,M, by ="row.names")

#testing only crans
marketclass <- candidates[(candidates$MarketClass == 'cran'),]
candGeno <- as.matrix(marketclass[,-c(1:3)]) #pull geno only

sampleIDs <- marketclass[,3]
GEBVs <- as.data.frame(candGeno %*% markerEffects)

GEBVs <- cbind(sampleIDs, GEBVs)
colnames(GEBVs) <- c("sampleIDs", "GEBVs")
Top <- GEBVs %>% arrange(desc((GEBVs)))
Top

write.csv(Top, "cranGEBVsJan2023")
