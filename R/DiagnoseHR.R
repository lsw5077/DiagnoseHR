#' Sensitivity diagnostics for home range estimates
#' 
#' hrDiag, assesses the role of individual relocations in driving estimate sensitivity 
#' by displaying how estimated home range size changes as individual relocations are iteratively removed.
#' 
#' @param locs is a SpatialPointsDataFrame containing coordinates(x,y) and 
#' a column of individual identifiers (ID). 
#' @param ID is a column from locs, and will be used for referencing individual-level diagnostics. 
#' @param type Defaults to 'mcp'. Options include 'mcp', 'kud' and 'LoCoH'
#' @param ... further arguments, including optional home range arguments from adeHabitatHR
#' @keywords home range
#' @export
#' @examples 
#' hrDiag(locs, type = 'mcp', ID = locs$ID)
#' 



hrDiag <- function(locs = NULL,
                   type = "mcp",
                   ID = locs$ID,
                   levels = 90,
                   percent = 95,
                   a = NULL,
                   unin = "m",
                   unout = "ha",
                   plotit = F,
                   grid = 60,
                   same4all = F,
                   kern = "bivnorm",
                   extent = 1,
                   duplicates = "random",
                   hlim = c(0.5, 1.5),
                   arange = NULL,
                   amount = "NULL",
                   h = "href"){
  

  
locs <- SpatialPointsDataFrame(locs, coords = locs[c("x","y")])


  
ID <- as.factor(ID)
locs$ID <- ID # assign loc ID as ID vector specified above. Will use for counting in loops.

base <- data.frame(ID = levels(ID), base_HR = NA) # make df to hold base hr results
result <- data.frame(ID = NA, loc_number = NA, HR_size = NA) # make result df to hold sample HR esitmats


  for (i in levels(ID)){ #our input df MUST have a column called "ID"
    
    subi <- locs[locs$ID == i,]

    base[base$ID == i,2] <- ifelse(type == "kud",
                   as.numeric(kernel.area(kernelUD(SpatialPoints(subi)))[paste(levels)], 
                              h = h,
                              grid = grid,
                              same4all = same4all,
                              kern = kern,
                              extent = extent,
                              boundary = NULL,
                              unin = unin,
                              unout = unout), # hr
                   ifelse (type == "LoCoH",
                           (as.numeric(LoCoH.a.area(SpatialPoints(subi),
                                                    arange = arange,
                                                    percent = percent,
                                                    unin = unin,
                                                    unout = unout,
                                                    duplicates = duplicates,
                                                    amount = amount))),
                           as.numeric(mcp.area(SpatialPoints(subi), 
                                               percent = percent,
                                               unin = unin,
                                               unout = unout,
                                               plotit = plotit)))) # set mcp as default
    
 
    
    for (j in 1:nrow(subi)){
            
      HR_size <- ifelse(type == "kud",
             as.numeric(kernel.area(kernelUD(SpatialPoints(subi[-j,])))[paste(levels)],
                        h = h,
                        grid = grid,
                        same4all = same4all,
                        kern = kern,
                        extent = extent,
                        boundary = NULL,
                        unin = unin,
                        unout = unout), # hr
             ifelse (type == "LoCoH", 
                     (as.numeric(LoCoH.a.area(SpatialPoints(subi[-j,]), 
                                              arange = arange,
                                              percent = percent,
                                              unin = unin,
                                              unout = unout,
                                              duplicates = duplicates,
                                              amount = amount))),
                     as.numeric(mcp.area(SpatialPoints(subi[-j,]),
                                         percent = percent,
                                         unin = unin,
                                         unout = unout,
                                         plotit = plotit))))    
      
     temp <- data.frame(ID = i, loc_number = j, HR_size = HR_size)
      
     result <- rbind(result, temp)
      
      }
      
  }
    # 
    result <- result[-1,]
    
    base$base_HR <- as.numeric(as.character(base$base_HR))
    base$ID <- as.factor(base$ID)
    result$ID <- as.factor(result$ID)
    
# Ok, maybe arrange by HR_size, assign index based on HR size
# 75% of the data are less than or equal to
# 25% of the data are less than or equal to
    
    
    result %>% as.data.frame() %>%
      left_join(base, by = "ID") %>% 
      mutate(deviation = (HR_size - base_HR)/base_HR) %>% as.data.frame() %>%
      group_by(ID) %>% 
      mutate(leverage = deviation/max(abs(deviation))) -> result
    

# result here should be a tibble, already grouped df      
      
      result %>%
      arrange(ID, HR_size) %>%
      mutate(mean_hr = mean(HR_size, na.rm = T),
      se = sd(HR_size)/length(HR_size),
      min = quantile(HR_size, method = 8)[[1]],
      q25 = quantile(HR_size, method = 8)[[2]],
      q50 = quantile(HR_size, method = 8)[[3]],
      max = quantile(HR_size, method = 8)[[4]],
      N = max(loc_number),
      hr_method = type) %>% 
      as.data.frame() %>% 
      mutate(method = type) %>%
      dplyr::select(ID, N, base_HR, mean_hr, se, min, q25, q50, max, hr_method) %>% distinct() -> result_tab
      
#sensitivity.plots <- list()  # create an empty list to save sensitivity plots
#leverage.plots <- list()   # create an empty list to save leverage plots in
  
for (i in levels(ID)){
  
  plot_result <- result[result$ID == i,]
   
  plot(plot_result$loc_number, plot_result$HR_size, pch = 19,
            main = paste("Home range sizes of individual", i, sep = " "),
            xlab = "Location number removed", ylab = "Home range size")
  lines(plot_result$loc_number, plot_result$HR_size)

  #sensitivity.plots[[i]] <- recordPlot() # save plot to sensitivity.plots list 
  #invisible(dev.off())
  
  hist(plot_result$leverage,
  main = paste("Leverage distribution of individual", as.character(i), sep = " "),
  xlab = "Leverage",
  ylab = "Frequency",
  breaks = 10)
  
#  leverage.plots[[i]] <- recordPlot() #save plot to leverage.plots list
# invisible(dev.off())
  
}
    
  out <<- list(result = result, result_tab = result_tab 
               #, sensitivity.plots = sensitivity.plots, 
               #leverage.plots = leverage.plots
              )
  print(out)
        
}

 
#' Asymptote diagnostics for home range estimates
#' 
#' The HRasym function allows users to examine whether and at what point the home range estimate
#' reaches an asymptote by iteratively adding relocations to a base subsample. 
#' The graph created by the HRasym function depicts the home range estimate for each iteration.
#' 
#' @param locs is a SpatialPointsDataFrame containing coordinates(x,y) and 
#' a column of individual identifiers (ID). 
#' @param ID is a column from locs, and will be used for referencing individual-level diagnostics. 
#' @param type Defaults to 'mcp'. Options include 'mcp', 'kud' and 'LoCoH'
#' @param ... further arguments, including optional home range arguments from adeHabitatHR
#' @keywords home range
#' @export
#' @examples 
#' hrDiag(locs, type = 'mcp', ID = locs$ID)
#' 

  
hrAsym <- function(locs = NULL,
                   type = "mcp",
                   ID = NULL,
                   levels = 90,
                   percent = 95,
                   a = NULL,
                   unin = "m",
                   unout = "ha",
                   plotit = F,
                   grid = 60,
                   same4all = F,
                   kern = "bivnorm",
                   extent = 1,
                   duplicates = "random",
                   hlim = c(0.5, 1.5),
                   arange = NULL,
                   amount = "NULL",
                   h = "href"){
  

 
  # get data
  
  ID <- as.factor(locs$ID)
  locs$ID <- ID
  locs <- SpatialPointsDataFrame(locs, coords = locs[c("x","y")])
  
  # make result df
  
  min_locs <- ifelse (type == "mcp", 5, 9)
  
  result <- as.data.frame(locs$ID) %>%
            mutate(ID = locs$ID) %>% dplyr::select(ID)%>%
            group_by(ID) %>%
            mutate(n_locs = row_number()) %>%
            filter(n_locs > min_locs) %>% mutate(HR_size = NA)
  
  asymptote.plots <- list()   # create an empty list to save asymptote plots to  
  
  # assign results in loop
  
  for (i in levels(ID)){
    
    plotHR <- locs[locs$ID == i,]
       
    for (j in 10:length(plotHR)){
      
      render <- plotHR[1:j,]
      
      result$HR_size[result$ID == i & result$n_locs == j] <- ifelse(type == "kud",
               as.numeric(kernel.area(kernelUD(SpatialPoints(render)))[paste(levels)], 
                          h = h,
                          grid = grid,
                          same4all = same4all,
                          kern = kern,
                          extent = extent,
                          boundary = NULL,
                          unin = unin,
                          unout = unout), # hr
               ifelse (type == "LoCoH",
                       (as.numeric(LoCoH.a.area(SpatialPoints(render),
                                                arange = arange,
                                                percent = percent,
                                                unin = unin,
                                                unout = unout,
                                                duplicates = duplicates,
                                                amount = amount))),
                       as.numeric(mcp.area(SpatialPoints(render), 
                                           percent = percent,
                                           unin = unin,
                                           unout = unout,
                                           plotit = plotit)))) # set mcp as default  
      
    }
    
    
    plot(result$n_locs[result$ID == i], result$HR_size[result$ID == i], pch = 19,
         main = paste("Home range size of individual", i, sep = " "),
         xlab = "Number of locations included in home range estimate", ylab = "Home range size")
    lines(result$n_locs[result$ID == i], result$HR_size[result$ID == i])

    #asymptote.plots[[i]] <- recordPlot() # save plot to asymptote.plots list
    #invisible(dev.off())
    
    }
  
  out <<- list(result = result) #, asymptote.plots = asymptote.plots
  print(out)
  
}
