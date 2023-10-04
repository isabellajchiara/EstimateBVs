load required packages
```{R}
library(caret)
library(ranger)
library(tidyverse)
library(BMTME)
library(e1071)
library(randomForest)
```

prepare data for model building
```{R}
data <- as.data.frame(read.csv("SNPdata.csv"))
data <- data %>%  select(where(~ n_distinct(.) > 1))
data <- na.omit(data)

M <- data[,-c(1,2)]
G <- A.mat(M) #rrblup function to calculate GRM

y <- data[,2]

data <- cbind(y,G)

colnames(data) <- paste("ID",1:(ncol(y) + ncol(G)), sep="")

```



create cross validation strategy 

```{R}
control <- trainControl(method='repeatedcv', 
                        number=10, ##will test 10 different values for mtry (number of variables for splitting) ##
                        repeats=3,
                        search = "random")  
```

build model

```{R}
rf_fit = train(ID1 ~ ., 
               data = data, 
               method = "rf",
               tuneLength = 10,
               trControl=control) ## search a random tuning grid 
```


predict GEBVs

```{R}
candidates <- read.csv("blackbeans.csv")
candidates <- as.matrix(candidates[,-1])
M <- A.mat(candidates)

GEBVs <- as.numeric(predict(rf_fit, candidates))
GEBVs <- as.data.frame(GEBVs)

```

Apply Sample IDs to GEBVs and order by value

```{R}
sampleIDs <- data[,1]

GEBVs <- cbind(sampleIDs, GEBVs)
colnames(GEBVs) <- c("IDs", "GEBVs")
Top <- GEBVs %>% arrange(desc((GEBVs)))

Top[1:5,]
```
