# DiagnoseHR: Sensitivity diagnostics for home range estimates 

Diagnose HR is an R package that provides functions for assessing sensitivity of home range estimates calculated using minimum convex polygon (MCP), local convex hull (LCH), and kernel utilization density (KUD) methods from the package [adehabitatHR](https://cran.r-project.org/web/packages/adehabitatHR/index.html). 

## Why assess home range sensitivity?

Variation in detection probability, sampling methods, and choice of home range analytical methods can influence the outcome of a home range analysis by influencing which and how many locations are included in the analysis. This can be a problem when the end goal of the analysis is to learn something about an ecosystem or plan a conservation action. When working on a home range analysis, however, it can be difficult to tell how sensitive a given home range estimate is to the sample size or the addition or removal of individual relocations. The functions below can help illustrate how sample size and the identities of relocations included in a home range analysis affect the outcome of that analysis. 

## Installing DiagnoseHR

DiagnoseHR can be installed using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) by copying and pasting the code below into an R script:

```{r}
install.packages("devtools") 
library(devtools) # install and load devtools

install_github("lsw5077/DiagnoseHR")
library(DiagnoseHR) # install DiagnoseHR from this github page
```

## Examining the importance of individual relocations

To see how much each relocation in your dataset is, you can use the function ```hrDiag()``` to recalculate your home range estimate leaving out one relocation at a time:

```{r}
load(sim_locs) # load the sample data that comes with the package

hrDiag(sim_locs) 

```







