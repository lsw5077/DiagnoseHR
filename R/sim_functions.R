
#square world functions:

#' Detection-explicit world building
#' 
#' make_world() builds a world according to your size and detection specifications

#' @param rows is the number of rows in the simulated world
#' @param columns is the number of columns in the simulated world
#' @param det_dist is the statistical distribution from which cell-specific
#' detection probabilities will be drawn. Options include a random uniform,
#' 'random', and beta distribution, 'beta'. Defaults to a random uniform between 0 and 1.
#' @param min_det is the minimum detection probability of any cell in the simulation
#' @param max_det is the maximum detection probability of any cell in the simulation
#' @param shape1 is the first shape parameter of the beta distribution.
#' @param shape2 is the second parameter of the beta distribution.
#' @keywords detection probability
#' @export
#' @examples
#' make_world(rows = 10, columns = 10, det_dist = 'random', min_det = 0, max_det = 1)

make_world <- function(rows,
                       columns,
                       det_dist = "random",# default distribution is the random uniform. Options = random uniform & beta  
                       min_det = 0,
                       max_det = 1,
                       shape1 = 10, # default shape makes symmetrical beta that approximates gaussian with mean = 0.5
                       shape2 = 10){
  
  expand.rows<- 1: rows  # expand the number of rows from 1 to number specified
  expand.cols<- 1: columns #expand the number of columns
  
  world<-expand.grid(y = expand.rows, x = expand.cols)  # make all combinations
  
  world$cell_id<- paste("id_", 1:nrow(world), sep="") # add a column of cell ids
  world$cell_prob<- NA # Default a column of probs with NAs
  
  if (det_dist != "random" & det_dist != "beta"){
    stop('Please choose a valid detection distribution') # error handling for invalid detection distributions
  }
  else if (det_dist == "random") {
    world$cell_prob <- runif(nrow(world), min_det, max_det) # if user selects random uniform, assign random uniform detection values
  }  
  else if (det_dist == "beta") {
    world$cell_prob <- rbeta(nrow(world), shape1, shape2) # if user selects beta, assign random beta detection values
  }  
  
  
  if (min_det < 0) {
    stop("Minimum detection probability must be greater than or equal to 0") # keep detection probs between 0 & 1
  } 
  else if (max_det > 1) {
    stop ("Maximum detection probability must be less than or equal to 1")
  }
  else {
    return(world)
  }
}


#' Populating simulated world with simulated organisims
#' 
#' The populate_world function initiates simulated organisms in a DiagnoseHR world
#' 
#' @param num_organisms is the number of simulated organisms to init
#' @param worldToPopulate is a make_world() detection world
#' @param num_per_cell is the maximum number of oganisms which may share a cell
#' @param replace cells may be sampled with or without replacement
#' @keywords correlated random walk
#' @export
#' @examples
#' populate_world(num_organisms = 3, worldToPopulate = world, num_per_cell = 1, replace = F)

populate_world <- function(num_organisms, # desired number of organisms
                           worldToPopulate, # the world to put them in
                           num_per_cell, # max number of orgs that can occupy a cell
                           replace = FALSE){ # no sampling with replacement

  potential_orgs <- worldToPopulate[rep(1:nrow(worldToPopulate), each=num_per_cell),] # gives you world with number of open slots = world * orgs per cell. 
  # If 4 orgs can fit on on a cell, then num_per_cell should be 4
  
  org_locations <- potential_orgs[sample(1:nrow(potential_orgs), num_organisms, replace=replace),] # define organisms by drawing locations from the potentials
  
  org_locations$org <- 1 # assign each org a population value, each one counts for 1
  
  loc_totals <- org_locations %>%
    group_by(x,y) %>%
    summarise(N=sum(org)) # count total number of orgs per cell
  
  square_world2 <- full_join(worldToPopulate, loc_totals, by=c("x", "y")) # returns starting world with organism counts
  
  square_world2$N[is.na(square_world2$N)] <- 0 # replace NA's w/0's
  
  return(square_world2)
}



# ok, we're going to be a bit creative with the homerange size argument. Users are going to expect
# a homerange in square units, but the functions need a radius. So we're going to convert from hr size using a = pi*r^2

run.walk<- function(Nsteps, # number of steps
                    x1, 
                    y1, 
                    homerange.radius, 
                    mu, 
                    rho,
                    wei.shape = 2,
                    wei.scale = homerange.radius,
                    site.fidelity){

  steplength <- rweibull(Nsteps, wei.shape, wei.scale) # draw steplengths from weibull dist
  
  thetaz <- suppressWarnings(rwrappedcauchy(Nsteps, mu = mu, rho = rho)) # draw directions
  uniformz <- runif(Nsteps, 0,1) # step-specific 
  
  walk.valz <- data.frame(step = 1:Nsteps,
                         steplength = steplength,
                         thetaz = thetaz, 
                         detect.prob = uniformz, 
                         x = NA, 
                         y = NA)
  
  walk.valz <- rbind(data.frame(step = 0,
                                steplength = 0, 
                                thetaz = 0, 
                                detect.prob = runif(1, 0,1),
                                x = x1, 
                                y = y1), 
                     walk.valz) # walk.valz gives you a movement trajectory for simulated orgs
  
  walk.valz <- find.step(walk.valz, 
                         homerange.radius,
                         site.fidelity) #but seems to need find.step, and homerange.size, which we define below
  
  return(walk.valz)
}


# distance finding funciton, helper, not for user

find.distance <- function(x1,
                          y1,
                          x2,
                          y2){
  distance<- sqrt((x1 - x2)^2 + (y1 - y2)^2) # simply finds distance between cells
  return(distance)
}


# direction finding function, helper, not for user

find.step<- function(walk.valz,
                     homerange.radius,
                     site.fidelity){
  
  walk.valz$new.dist<- NA
  walk.valz$prop.dist <- NA
  
  for(i in 2:nrow(walk.valz)){
    
    dX <- walk.valz$steplength[i]*cos(walk.valz$thetaz[i]) # calculate x distance based on trajectory
    dY <- walk.valz$steplength[i]*sin(walk.valz$thetaz[i]) # calculate y distnce based on trajectory
    
    potential.step.x <- as.numeric(walk.valz$x[i - 1] + dX) # land on potential x vector
    potential.step.y <- as.numeric(walk.valz$y[i - 1] + dY) # land on potential y vector
    
    walk.valz$new.dist[i] <- find.distance(x1 = walk.valz$x[1],
                                           y1 = walk.valz$y[1], 
                                           potential.step.x, 
                                           potential.step.y) # calculate distance between current loc and potential step
    
    walk.valz$prop.dist[i] <- walk.valz$new.dist[i]/homerange.radius # calculate how much of the home range radius the step would be
    
    probz <- site.fidelity*walk.valz$prop.dist[i]*2/(1+ site.fidelity*walk.valz$prop.dist[i])
    
    if(probz > runif(1,0,1)){
      #  If probability exceeds a random uniform, organism turns back toward origin point
      walk.valz$thetaz[i]<- atan2((walk.valz$y[1]-potential.step.y), (walk.valz$x[1]-potential.step.x))
    }
    
    
    walk.valz$x[i]<-as.numeric(walk.valz$x[i - 1] + walk.valz$steplength[i] * cos(walk.valz$thetaz[i]))
    walk.valz$y[i]<-as.numeric(walk.valz$y[i - 1] + walk.valz$steplength[i] * sin(walk.valz$thetaz[i]))
  }
  
  return(walk.valz[,c("step","detect.prob","x","y")])
}



#' Simulate correlated random walks using the simulated population (pop_world) and the detection space(myworld)
#' 
#' 
#' @param pop_world is the simulated population, created by populate_world()
#' @param myworld is a detection world, created by make_world()
#' @param world.type worlds may be 'open' such that orgs can leave the sample space, or 'closed' such that they cannot
#' @param Nsteps the number of timesteps in which the organisms can make movements
#' @param homerange.type Home ranges may be either 'fixed' with a defined home range for all organisms, or 'random' such that each individual draws a home range size in each iteration
#' @param homerange.size Home range size in square units of individuals with a fixed home range
#' @param mu the mu parameter of the wrapped cauchy distribution determining step direction 
#' @param rho the rho parameter of the wrapped cauchy distribution determining step direction 
#' @param wei.shape the shape parameter of the weibull distribution determining step length
#' @param wei.scale the scale parameter of the weibull distribution determing step length
#' @param site.fidelity the degree of site fidelity (0-1)
#' @keywords correlated random walk
#' @export
#' @examples
#' move_critters(pop_world = pop_world, myworld = world, world.type = 'closed', Nsteps = 10, homerange.size = 10, mu = 0, rho = 0, site.fidelity = 1)


move_critters <- function(pop_world, 
                          myworld, 
                          world.type="closed",
                          Nsteps, 
                          homerange.type="fixed", 
                          homerange.size, 
                          mu = 0, 
                          rho = 0,
                          wei.shape = 2,
                          wei.scale = sqrt(homerange.size/pi),
                          site.fidelity = 1){
  
  if (site.fidelity > 1 | site.fidelity < 0) {
    stop 
    print("Please choose a site fidelity between 0 and 1")}
  
  homerange.radius <- sqrt(homerange.size/pi)

  pop_world.red<-pop_world[pop_world$N>0,] #reduced population to cells with organisms
  
  pop_world.red<-pop_world.red[rep(1:nrow(pop_world.red),pop_world.red$N),c("x","y","cell_id")] # expands the cells for each organism
  
  rownames(pop_world.red)<-NULL
  
  num.org<-nrow(pop_world.red)  #calculate the number of unique organisms
  
  org_stor<-  as.data.frame(matrix(NA,Nsteps*num.org,5))  # create a data.frame to stor organism movement
  names(org_stor)<-c("step","detect.prob","x","y","org_id") # rename column headers
  
  for(i in 1:num.org){
    if(homerange.type=="fixed"){
      run_steps<-run.walk(Nsteps, 
                          x1=pop_world.red$x[i], 
                          y1=pop_world.red$y[i], 
                          homerange.radius = homerange.radius, # homerange fixed for all organisms
                          mu = mu,
                          rho=rho,
                          site.fidelity = site.fidelity)  
      
    }
    if(homerange.type=="random"){
      run_steps<-run.walk(Nsteps, 
                          x1=pop_world.red$x[i],
                          y1=pop_world.red$y[i], 
                          homerange.size=rgamma(1,homerange.radius), # mean homerange size is equal to homerange size (poisson distribution)
                          mu = mu, 
                          rho=rho,
                          site.fidelity = site.fidelity) 
      
    }
    
    run_steps$org_id<-paste("org",i,sep="_")
    
    org_stor[(Nsteps*(i-1)+i):((Nsteps*i +i)),] <- run_steps
    
  }
  
  if(world.type=="closed"){
    min.y<- min(myworld$y)   #  Organism stays on edge, can not leave myworld
    max.y<- max(myworld$y)
    min.x<- min(myworld$x)
    max.x<- max(myworld$x)
    
    org_stor$y[org_stor$y<min.y]<-min.y
    org_stor$y[org_stor$y>max.y]<-max.y
    org_stor$x[org_stor$x<min.x]<-min.x
    org_stor$x[org_stor$x>max.x]<-max.x
    
  }
  
  org_stor$cell_x<-floor(org_stor$x) # round down to find cell number in x and y
  org_stor$cell_y<-floor(org_stor$y)
  
  org_stor<-left_join(org_stor,myworld, by=c("cell_x"="x", "cell_y"="y"))  # join to myworld to get cell specific detection
  
  if(world.type=="open"){  
    org_stor$cell_prob[is.na(org_stor$cell_prob)]<-0 # sets cell specific detection to 0 if organism moves off myworld
  }
  
  org_stor <- org_stor %>% 
    dplyr::mutate(do.detect = ifelse(detect.prob < cell_prob,1,0),
           step = step + 1)# create a binary value for detection. 1 = detect
  
  return(org_stor)
  
}



#' Sample a movement record from a correlated random walk in a detection landscape
#' 
#' 
#' @param world The detection world to sample
#' @param walk The correlated random walk to sample
#' @param n.cells The number of cells to sample per timestep
#' @param sample.steps The number of time steps to sample
#' @param replace.step Whether to sample with replacement within a time step
#' @param replace.world Whether to sample with replacement within the entire simulation
#' @param seed An optional seed
#' @keywords sampling, detection
#' @export
#' @examples
#' sample_world(world = world, walk = walk, sample.steps = seq(1,10,1), replace.step = F, replace.world = T, seed = NULL)

sample_world <- function(world, # the movement record to sample
                         walk, # the movement record
                         n.cells, # the total number of cells to be sampled per day
                         sample.steps, # a list of the number of timesteps to be sampled. 
                         replace.step = F, # replacement within a timestep
                         replace.world = T, #replace any cells in the whole sim
                         seed=NULL){
  
  sample_stor <- data.frame(cell_id = matrix(NA,n.cells*length(sample.steps)))
  
  set.seed(seed)
  
  for (i in 1:length(sample.steps)){
    
    if(replace.world == T) {
      available.cells <- world$cell_id # create object of available points. In this case is all the rows that we have initially
    }
    
    if(replace.world == F){
      available.cells <- world$cell_id[-world$cell_id %in% sample_stor$cell_id] # world without replacement
    }
    
    sample_stor$cell_id[seq(((i-1)*(n.cells)) + 1, i*n.cells,1)] <- sample(available.cells, n.cells, replace = replace.step) # sample the available cells
    sample_stor$step[seq(((i-1)*(n.cells)) + 1, i*n.cells,1)] <- sample.steps[i] 
    
  }
  
  sample_walk <- left_join(sample_stor, walk, by = c("cell_id","step")) %>%
                 dplyr::select(-x,-y)
  
  detect_walk <- sample_walk[sample_walk$do.detect == 1,]
  org_sample <- list(cells_sampled = sample_stor, orgs_in_sample = na.omit(sample_walk), orgs_detected = na.omit(detect_walk))
  
  return(org_sample)
}   

