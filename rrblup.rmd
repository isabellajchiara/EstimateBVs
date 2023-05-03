load required packages
```{R}
library(rrBLUP)
library(BMTME)
library(tidyverse)
library(readxl)
```

load data and prepare matrices

data <- as.data.frame(read_xlsx("trainingData.xlsx", sheet= 2))

```{R}
geno <- data[,-c(1:15)]
pheno <- as.data.frame(data[,6])
pheno <- as.matrix(pheno[1:269,])
ids <- as.data.frame(data[,1])

geno[geno=="--"] <- NA
genoVar <- geno %>%  select(where(~ n_distinct(.) > 1))
genoVar <- genoVar[1:269,]
M <- as.matrix(genoVar %>%select_if(~!any(is.na(.))))


```

Build model and predict marker effects
```{R}

EBVans <-mixed.solve(y=pheno, Z=M, K=NULL, SE=FALSE, return.Hinv=FALSE)

markerEffects <- EBVans$u

```


candidates <- read_xlsx("Candidates.xlsx", sheet=2)
```{R}
marketclass <- candidates[(candidates$MarketClass == 'cran'),]
candGeno <- as.matrix(marketclass[,-c(1,2)])
sampleIDs <- marketclass[,2]
```

```{R}
GEBVs <- as.data.frame(candGeno %*% markerEffects)
```

Apply sample IDs to GEBVs and order by value
```{R}

GEBVs <- cbind(sampleIDs, GEBVs)
colnames(GEBVs) <- c("IDs", "GEBVs")
Top <- GEBVs %>% arrange(desc((GEBVs)))

cranParentsReclustSNP <- as.data.frame(Top[1:5,1])
```
