
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
#' @examples
#' populate_world(num_organisms = 3, worldToPopulate = world, num_per_cell = 1, replace = F)

populate_world <- function(num_organisms, # desired number of organisms
                           worldToPopulate, # the world to put them in
                           num_per_cell, # max number of orgs that can occupy a cell
                           replace = FALSE){ # no sampling with replacement
  require(dplyr) #### relocate
  
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



# 

run.walk<- function(Nsteps, # number of steps
                    x1, 
                    y1, 
                    wei_shape = 2,
                    wei_scale = 350,
                    homerange.size, 
                    mu, 
                    rho){
  require(circular) ## relocate
  
  steplength <- rweibull(Nsteps, wei_shape, wei_scale) # draw steplengths from weibull dist
  
  thetaz <- suppressWarnings(rwrappedcauchy(Nsteps, mu = mu, rho = rho)) # draw directions
  uniformz <- runif(Nsteps, 0,1) # step-specific 
  
  walk.valz<- data.frame(step = 1:Nsteps,
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
  
  walk.valz <- find.step(walk.valz, homerange.size) #but seems to need find.step, and homerange.size, which we define below
  
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
                     homerange.size){
  
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
    
    walk.valz$prop.dist[i] <- walk.valz$new.dist[i]/homerange.size # calculate how much of the home range the step would be
    
    probz <- (0.9999 / (1 + exp(4.741 + -9.407*walk.valz$prop.dist[i]))) 
    
    if(probz > runif(1,0,1)){
      #  If probability exceeds a random uniform, organism turns back toward origin point
      walk.valz$thetaz[i]<- atan2((walk.valz$y[1]-potential.step.y), (walk.valz$x[1]-potential.step.x))
    }
    
    
    walk.valz$x[i]<-as.numeric(walk.valz$x[i - 1] + walk.valz$steplength[i] * cos(walk.valz$thetaz[i]))
    walk.valz$y[i]<-as.numeric(walk.valz$y[i - 1] + walk.valz$steplength[i] * sin(walk.valz$thetaz[i]))
  }
  
  return(walk.valz[,c("step","detect.prob","x","y")])
}


# Move the animals, record their movements

move_critters <- function(pop_world, 
                          myworld, 
                          world.type="closed",
                          Nsteps, 
                          homerange.type="fixed", 
                          homerange.size, 
                          mu = 0, 
                          rho = 0){
  require(dplyr)
  require(circular)
  
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
                          homerange.size=homerange.size, # homerange fixed for all organisms
                          mu = mu, rho=rho)  
      
    }
    if(homerange.type=="random"){
      run_steps<-run.walk(Nsteps, 
                          x1=pop_world.red$x[i],
                          y1=pop_world.red$y[i], 
                          homerange.size=rpois(1,homerange.size), # mean homerange size is equal to homerange size (poisson distribution)
                          mu = mu, 
                          rho=rho) 
      
    }
    
    run_steps$org_id<-paste("org",i,sep="_")
    
    org_stor[(Nsteps*(i-1)+i):((Nsteps*i +i)),]<-run_steps
    
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
  
  org_stor<-org_stor %>% 
    mutate(do.detect = ifelse(detect.prob<cell_prob,1,0))  # create a binary value for detection. 1 = detect
  
  return(org_stor)
  
}



# ok, so we have N_samples. Need to work in timing of samples too. 
# just need to add a time component to output sample grid
# join walk df to sample grid by x, y, time



# make simpler sampling cells with replacement


sample_myworld <- function(world, # the movement record to sample
                           #walk, # the movement record
                           n.grids, # the total number of grids sampled
                           size.x, # x dimensions of sampling grids
                           size.y, # y dimensions of sampling grids
                           sample.days, # a list of the days to be sampled
                           max.iter = 1000,
                           seed=NULL){
  
  #creates a function to generate random sampling grids in myworld 
  if(!is.null(seed)){
    set.seed(seed) # if seed is not null, set the seed to seed
  }
  
  s = 1 # counter for while loop
  t = 1 # counter for total iterations

  available.cells <- world$cell_id # create object of available points. In this case is all the rows that we have initially
  
  sample_stor <- data.frame(matrix(NA,(size.x*size.y*n.grids*length(sample.days)), 5)) #creates storage for sample cells, expanding the matrix by sample days
  names(sample_stor) <- c(names(world), "grid.id") # rename columns, make a step variable so we can match it up later
  
  
  while(s <= (n.grids*length(sample.days))){
    start.point <- world[sample(1:nrow(world), size=1, replace = T),] # sample an initial row
    
    if(start.point$cell_id %in% available.cells){ 
      sample.dir <- data.frame(UD = sample(c('U','D'), size=1), LR = sample(c('L','R'), size=1)) # sample directions, left, right, up, down
      
      if(sample.dir$UD =='U'& sample.dir$LR =='R'){
        poss.x<-seq(start.point$x, start.point$x + (size.x-1), by = 1) # if up and right, create vectors of positive x
        poss.y<-seq(start.point$y, start.point$y + (size.y-1), by = 1) # if up and right, create vectors of positive y
      }
      
      if(sample.dir$UD=='U'& sample.dir$LR =='L'){
        poss.x<-seq(start.point$x, start.point$x - (size.x-1), by = -1) # if up and right, create vectors of negative x
        poss.y<-seq(start.point$y, start.point$y + (size.y-1), by = 1) # if up and right, create vectors of positive y
      }
      
      if(sample.dir$UD=='D'& sample.dir$LR=='R'){
        poss.x<-seq(start.point$x, start.point$x + (size.x-1), by=1) # if up and right, create vectors of positive x
        poss.y<-seq(start.point$y, start.point$y - (size.y-1), by=-1) # if up and right, create vectors of negative y
      }
      
      if(sample.dir$UD=='D'& sample.dir$LR=='L'){
        poss.x<-seq(start.point$x, start.point$x - (size.x-1), by=-1) # if up and right, create vectors of negative x
        poss.y<-seq(start.point$y, start.point$y - (size.y-1), by=-1) # if up and right, create vectors of negative y
      }
      
      poss.cells<-expand.grid(x=poss.x, y=poss.y) # create sample grid of possible x-y combos
      poss.cells<-left_join(poss.cells, world, by=c("x","y")) # grab the part of the world that we intend to sample
      
      if(all(poss.cells$cell_id %in% available.cells)){
        sample.cells <- poss.cells
        sample.cells$sample_id< - s
        
        sample_stor[(((s-1)*(size.x*size.y))+1): (s*size.x*size.y),] <- sample.cells
        
        # start ((s-1)*(size.x*size.y))+1
        # end s*size.x*size.y
        
        available.cells <- available.cells[-which(available.cells %in% sample.cells$cell_id)] # remove sampled cells from available
        
        s <- s+1 # 
        
       } 
    }    
    t = t + 1
    if(t == max.iter){
      stop("You have reach the maximun number of available iterations")
    }
    
  }
  
  
  sample_stor$sample.id <- paste(sort(rep(seq(1, n.grids, 1), length(sample.days))),
                                 sort(rep(seq(1, size.x*size.y,1), length(sample.days)*n.grids)), sep = ".")
  sample_stor$step <- rep(sample.days, n.grids) # assign indexes 

  sample_walk <- left_join(sample_stor, walk, by = c("step", "x", "y", "cell_id", "cell_prob")) 
  detect_walk <- sample_walk[sample_walk$do.detect == 1,]
  spacetime_sample <- list(Grids_sampled = sample_stor, orgs_in_sample = na.omit(sample_walk), orgs_detected = na.omit(detect_walk))
  
  
  return(spacetime_sample)
  
  
  #sample_walk <- left_join(sample_stor, walk, by = c("x", "y", "cell_id", "cell_prob", "step"))
  #walk_sample <- list()
  
  
}


