## Introduction

A R package for biological experiments analysis. Beta version now.

## Installation

```R
# install luckyExperiment
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('huangwb8/luckyExperiment')
```

library `luckyExperiment` package, just do: 

```R
library(luckyExperiment)
```

Thus, every functions next in `luckyExperiment` can only be available.

## Functions

### 1. FastQPCR

#### package

install other needed packages

```R
package.need <- c('tidyr','readxl')
for( i in package.need){
  if (!requireNamespace(i, quietly = TRUE)){
    install.packages(i)
    library(i,character.only = T)
  } else {
    library(i,character.only = T)
 }
} 
```

or you can just simply use `luckyExperiment::Plus.library` to do this:

```R
Plus.library(c('tidyr','readxl'))
```

#### example data

```R
# a RT-PCR data for analysis
data <- system.file("extdata","testData_qPCR.xlsx", package = "luckyExperiment") %>% read_xlsx %>% as.data.frame
View(data)
```

#### FastQPCR analysis

```{r}
result_qpcr <- FastQPCR(data,
                        sample = "samples", 
                        marker = "markers",
                        bioRepeat = "biorepeat", 
                        parallelRepeat = "parepeat",
                        group = "groups", 
                        group.control = "sh-Negative", 
                        internal = "GAPDH",
                        value = "Ct", 
                        plot.type = 1, # one marker
                        palette = NULL,
                        size = 20, 
                        label.position  = c(1,1.5),
                        label.type = "p.format",
                        method =  "anova",
                        x.title = "Groups", 
                        y.title = "The Relative Expression of Genes",
                        legend.position = "top")

# result
View(result_qpcr$Data$metadata)
View(result_qpcr$Data$statistc$whole)
View(result_qpcr$Data$statistc$pair) 
```