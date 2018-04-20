# DiagnoseHR: Sampling simulation and sensitivity diagnostics for home range estimates (and beyond!)

Diagnose HR is an R package that provides functions for simulating detection-censored animal movement data and conducting home range sensitivity analysis. DiagnoseHR can facilitate sensitivity analysis of home ranges constructed using minimum convex polygon (MCP), local convex hull (LCH), and kernel utilization density (KUD) methods from the package [adehabitatHR](https://cran.r-project.org/web/packages/adehabitatHR/index.html). 

## Why consider detection probability?

Variation in detection probability, sampling methods, and choice of home range analytical methods can influence the outcome of a spatial analysis by influencing which and how many locations are included in the analysis. This can be a problem when the end goal of the analysis is to learn something about an ecosystem or plan a conservation action. When working on a movement or home range analysis, however, it can be difficult to tell how sensitive a given estimate is to the sample size or the addition or removal of individual relocations. The functions in DiagnoseHR can help simulate data that reflects realistic sampling and detection constraints. The home range sensitivity functions can help illustrate how sample size and the identities of relocations included in a home range analysis affect the outcome of that analysis. 

## Installing DiagnoseHR

DiagnoseHR can be installed using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) by copying and pasting the code below into an R script:

```{r}
install.packages("devtools") 
library(devtools) # install and load devtools

install_github("lsw5077/DiagnoseHR")
library(DiagnoseHR) # install DiagnoseHR from this github page
```
## Simulating data

Simulating data with DiagnoseHR occurs in four steps that reflect process by which animal movement data is created and collected: landscape formation, population establishment, animal movement, and sampling

### Landscape formation

We start by using the ```make_world()``` function to create a landscape that varies in detection probability, just like a real landscape. Our world, ```world```, is a matrix of dimensions 10x10. Its detection distribution is defined by the flexible beta distribution. By setting the shape1 and shape2 arguments to the same number, we approximate a normal distribution bounded between 0 and 1. When we graph ```world```, we can see that the tile colors reflect the probability of observing an organism in the space given its presence.

```{r}
world <- make_world(rows = 10, columns = 10, det_dist = "beta", shape1 = 10, shape2 = 10)

ggplot()+
      geom_tile(data = world, aes(x = x, y = y, fill = cell_prob))+
      coord_cartesian(expand = F)+
      theme_bw()

```

![Detection_world](/images/world_plot.png)

### Populating the world

Now that we have our simulation space, ```world```, we need an organism to observe. We populate the simulation space with one organism using the ```populate_world()``` function. When we look at our plot now, we can see our organism in the simulation space. 

```{r}

popworld <- populate_world(1, # spawn one organism
                           world, # in world
                           1, # with a maximum of one organism per cell
                           replace = T) # sampling with replacement

ggplot()+
  geom_tile(data = popworld, aes(x = x, y = y, fill = cell_prob))+
  geom_point(data = popworld[popworld$N > 0,], aes(x = x, y = y), size = 20, color= "black") +
  geom_point(data = popworld[popworld$N > 0,], aes(x = x, y = y), size = 10, color= "red") +
  coord_cartesian(expand = F)+
  theme_bw()
```

![Popworld](/images/popworld.png)

### Simulating organism movement

Our simulated organism will now take a walk around the simulation space, which we will direct using the ```move_critters()``` function.

```{r}
walk <- move_critters(popworld, # population to move
                      world, # simulation space to move through
                      world.type = "closed", # can the orgs move off the simulation space ("open") or not? ("closed")
                      Nsteps = 20, # Number of time steps during which the organism can move
                      homerange.size = 10, # the organism's home range size, if fixed
                      homerange.type = "fixed", # do organisms get the same home range size ("fixed") or draw each iteration ("random")?
                      mu = 0, # direction distribution, from a wrapped cauchy distribution
                      rho = 0)
                      
ggplot()+
  geom_tile(data = world, aes(x = x, y = y, fill = cell_prob))+
  geom_line(data = walk, aes(x = x, y = y), size = 2)+
  geom_point(data = walk, aes(x = x, y = y), size = 20, color= "black") +
  geom_point(data = walk, aes(x = x, y = y), size = 10, color= "red") +  coord_cartesian(expand = F)+
  theme_bw()                      
                      
```
![Correlated_random_walk](/images/walk.png)

### Sampling organism movement

So far, we have created a simulation space with varying detection, and a simulated organism that travels around the simulation space. Now let's look at what happens when impose some realistic sampling constraints. Using the ```sample_world()``` function, we can see what happens to our organism movement record when we sample half the cells every other day for the duration of the study period. When we compare our sampling plot to our plot of the organism's walk, we can see that, even though the path is similar, we missed some observations because we didn't sample the right time step or detection probability was too low. 

```{r}

sample <- sample_world(world, # The simulation space in which to sample
                       walk, # The organism movement record to sample
                       n.cells = 50, # the number of cells to sample per time step
                       sample.steps = seq(2, 20, 2), # the time steps to sample
                       replace.step = F, # don't allow replacement within timesteps
                       replace.world = T, # do allow replacement between timesteps,
                       seed = NULL) # We don't set a seed, although we could if we want consistent results
                       
 ggplot()+
  geom_tile(data = world, aes(x = x, y = y, fill = cell_prob))+
  geom_line(data = walk, aes(x = x, y = y), size = 2)+
  geom_point(data = sample$orgs_detected, aes(x = cell_x, y = cell_y), size = 20, color= "black") +
  geom_point(data = sample$orgs_detected, aes(x = cell_x, y = cell_y), size = 10, color= "red") +  coord_cartesian(expand = F)+
  theme_bw()                      
                       
```
![Sample](/images/sample.png)

## Assessing home range estimate sensitivity

### Examining the importance of individual relocations

hrDiagnose requires a relocation dataset with an individual identification column "ID", x coordinate "x", and y coordinate "y". To see how much each relocation in your dataset is, you can use the function ```hrDiag()``` to recalculate your home range estimate leaving out one relocation at a time:

```{r}
load(sim_locs) # load the sample data that comes with the package
hrDiag(locs = sim_locs, ID = sim_locs$ID) 
```
The resulting elasticity plot shows how the size of the home range changes as each relocation is removed. A home range estimate where each relocation has a large effect on the size of the estimate will produce an elasticity plot with large, jagged fluctuations, while a less sensitive estimate will appear smoother. HRdiag also calculates the leverage of each relocation, how much each point contributes to the difference between the estimate with all the relocations and each N-1 recalculated home range. A leverage histogram where each relocation contributes a small amount will have a large number of observations clustered around 0, while an estimate where a few relocations have a large effect will produce a leverage histogram with observations distributed away from 0.

![Elasticity plot](/images/9_desert.png) ![Leverage plot](/images/9_desert_lev.png)

### Assessing home range asymptotes

Home ranges are traditionally assumed to reach an asymptote as the number of observations increases. There are many reasons a home range may not reach an asymptote, and a home range that is not asympototic may still be useful. However, the interpretation and subsequent applicaiton of the estimate may change. The function ```hrAsym()``` iteratively adds relocations to an initial subsample to show if and when a home range reaches an asymptote:

```{r}
hrAsym(locs = sim_locs, ID = sim_locs$ID)
```
![Asymptote plot](/images/asym.png)

### Disclaimer

Although these data have been processed successfully, no warranty expressed or implied is made regarding the display or utility of the data on any other system or for general or scientific purposes, nor shall the act of distribution constitute any such warranty. The USGS or the U.S. Government shall not be held liable for improper or incorrect use of the data described and (or) contained herein.
