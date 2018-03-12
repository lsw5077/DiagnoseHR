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

hrDiagnose requires a relocation dataset with an individual identification column "ID", x coordinate "x", and y coordinate "y". To see how much each relocation in your dataset is, you can use the function ```hrDiag()``` to recalculate your home range estimate leaving out one relocation at a time:

```{r}
load(sim_locs) # load the sample data that comes with the package
hrDiag(sim_locs) 
```
The resulting elasticity plot shows how the size of the home range changes as each relocation is removed. A home range estimate where each relocation has a large effect on the size of the estimate will produce an elasticity plot with large, jagged fluctuations, while a less sensitive estimate will appear smoother. HRdiag also calculates the leverage of each relocation, how much each point contributes to the difference between the estimate with all the relocations and each N-1 recalculated home range. A leverage histogram where each relocation contributes a small amount will have a large number of observations clustered around 0, while an estimate where a few relocations have a large effect will produce a leverage histogram with observations distributed away from 0.

![Elasticity plot](/images/9_desert.png) ![Leverage plot](/images/9_desert_lev.png)

## Assessing home range asymptotes

Home ranges are traditionally assumed to reach an asymptote as the number of observations increases. There are many reasons a home range may not reach an asymptote, and a home range that is not asympototic may still be useful. However, the interpretation and subsequent applicaiton of the estimate may change. The function ```hrAsym()``` iteratively adds relocations to an initial subsample to show if and when a home range reaches an asymptote:

```{r}
hrAsym(sim_locs)
```
![Asymptote plot](/images/asym.png)










