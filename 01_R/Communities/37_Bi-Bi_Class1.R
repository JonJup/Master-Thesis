### --- Bi-Bi Class 1 --- ###  

#------------#
#1. SETUP ----
#------------#
    
    # sets wd to location of sourced script. Only works when sourced! 
    setwd(dirname(parent.frame(2)$ofile))
    # The file collected_functions.R contains all the functions necessary to run this script. 
    source("Collected_Functions.R")
    
#----------------#
#2. PARAMETER ----
#----------------#

    # number of species 
    n <- 9
    
    # number of sampling sites per gradient 
    sites <- 10
    
    # number of noise variables 
    nerr <- 1
    
    # env is the grid of x and y gradient along which species abundances will be simulated.     
    env <- expand.grid(1:100,1:100)
    
    # Now the species' response pararmeters are set:
    # maximal abundance 
    c <-    list(
                 "B" = list(
                        list(c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100)),
                        list(c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100))
                 )
            )
      
    # optima                   
    opt <-  list(
                "B" = list(
                        list(c(5,25),c(5,25),c(5,25),c(35,55),c(35,55),c(35,55),c(75,95),c(75,95),c(75,95)),
                        list(c(5,25),c(35,55),c(75,95),c(5,25),c(35,55),c(75,95),c(5,25),c(35,55),c(75,95))
                ) 
               )
               
    # tolerance
    tol <- list ("B" = list ( list(c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6)),
                               list(c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6))))
    
    # combine parameters into list 
    para <- list(c = c, opt = opt, tol = tol)


#------------------#
# 3. Simulation ----
#------------------# 

    # This function simultes the abundances. For more details see the collected_functions script 
    out <- sim.re.md(n.grad = 2, re.type = list(rep(4, n), rep(4, n)), parameter = para, grad = env) 
    
    # Sample the big 100x100 dataset    
    names(out) <- paste("sp", 1:ncol(out), sep = "")
    
    # Sample the big 100x100 dataset
    sample <- take.sample(out, n.x = sites, n.y = sites, grid = env)
    
    # Extract new abundance data
    out.samp <- sample[["Abundance"]]
    
    # Extract new sample location 
    out.grid <- sample[["Sample grid"]]
    
    # Add noise variables 
    envn <- make.noise(n.grad = nerr, n.cor = .0, grid = out.grid, h = 100)