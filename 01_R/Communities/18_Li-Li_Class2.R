### --- Li-Li Class 2 --- ###  

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
    nerr <- 5
    
    # env is the grid of x and y gradient along which species abundances will be simulated.     
    env <- expand.grid(1:100,1:100)
    
    # Now the species' response pararmeters are set: 
    # linear parameter 
    lc = list(c(80,100,120,0.8,1.1,1.2),
              c(80,100,120,0.8,1.1,1.2))

    # combine all parameters into one list 
    para <- list(lc = lc)

#------------------#
# 3. Simulation ----
#------------------# 
    
    # This function simultes the abundances. For more details see the collected_functions script     
    out <- sim.re.md(n.grad = 2, re.type = list(c(1,1,1,0,0,0,1,1,1),c(0,0,0,1,1,1,1,1,1)), parameter = para, grad = env)

    # Sample the big 100x100 dataset
    sample <- take.sample(out, n.x = sites, n.y = sites, grid = env)
    
    # Extract new abundance data
    out.samp <- sample[["Abundance"]]
    
    # Extract new sample location 
    out.grid <- sample[["Sample grid"]]
    
    # Add noise variables 
    envn <- make.noise(n.grad = nerr, n.cor = .0, grid = out.grid, h = 100)