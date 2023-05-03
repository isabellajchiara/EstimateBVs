load required packages
```{R}
library(e1071)
library(tidyverse)
library(BMTME)
library(rrBLUP)
```

Prepare data for model building
```{R}
data <- as.data.frame(read.csv("SNPdata.csv"))
data <- data %>%  select(where(~ n_distinct(.) > 1))
data <- na.omit(data)

M <- data[,-c(1,2)]
G <- as.data.frame(A.mat(M)) #rrblup function to calculate GRM

y <- as.data.frame(data[,2])

data <- cbind(y,G)


colnames(data) <- paste("ID",1:(ncol(y) + ncol(G)), sep="")

```

Fit model and predict GEBVs
```{R}
svm_fit = svm(ID1 ~ ., data = data, kernel = "radial", cost = 10, scale = FALSE)
```

predict GEBVs for selection candidates

```{R}
candidates <- read.csv("blackbeans.csv")
candidates <- as.matrix(candidates[,-1])
M <- as.data.frame(A.mat(candidates))
colnames(M) <- paste("ID",3:(ncol(M)+1),sep="")
GEBVs <- predict(svm_fit, M)

```

Apply sample IDs to GEBVs and order by value

```{R}
sampleIDs <- data[,1]

GEBVs <- cbind(sampleIDs, GEBVs)
colnames(GEBVs) <- c("IDs", "GEBVs")
Top <- GEBVs %>% arrange(desc((GEBVs)))

Top <- as.data.frame(Top)
Top[1:5,]

```
