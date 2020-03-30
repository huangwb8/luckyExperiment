## Introduction

A R package for biological experiments analysis. Beta version now.

## Installation

```R
# install luckyExperiment
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github('huangwb8/luckyExperiment')
```

## Functions

### 1. FastQPCR

#### example data

```R
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