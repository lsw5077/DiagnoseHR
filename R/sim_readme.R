library(tidyverse)
library(RColorBrewer)

world <- make_world(rows = 50, columns = 50, 
                    det_dist = "random", min_det = 0, max_det = 1)

ggplot()+
      geom_tile(data = world, aes(x = x, y = y, fill = cell_prob))+
      theme_bw()

popworld <- populate_world(3, world, 2, replace = T)

ggplot()+
  geom_tile(data = popworld, aes(x = x, y = y, fill = cell_prob))+
  geom_point(data = popworld[popworld$N > 0,], aes(x = x, y = y), size = 2, color= "red") +
  theme_bw()

walk <- move_critters(popworld,
                      world,
                      world.type = "closed",
                      Nsteps = 10,
                      homerange.size = 10,
                      homerange.type = "fixed",
                      mu = 0,
                      rho = 0)

ggplot()+
  geom_tile(data = world, aes(x = x, y = y, fill = cell_prob))+
  geom_line(data = walk, aes(x = x, y = y, color = org_id), size = 2)+
  geom_point(data = walk, aes(x = x, y = y, color = org_id),  size = 2) +
  theme_bw()


sample_myworld(world, n.grids = 4, 
               sample.days = sample.days,  
               size.x = 10, size.y = 10, max.iter = 200000)


sample.days <- seq(0, 3, 1)

test <- data.frame(seq(0, 100, 2))
test$sample <- sample.days

#same places each day, then 

