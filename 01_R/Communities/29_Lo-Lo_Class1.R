### --- Lo-Lo Class 1 --- ###  

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
    # maximal abundances 
    c <-    list("L" = list(rep(100, n), rep(100, n)))
    
    # Optima
    opt <-  list("L" = list( 
        c(20, 50, 80, 20, 50, 80, 20, 50, 80),
        c(20, 20, 20, 50, 50, 50, 80, 80 ,80)))
    
    # logistic parameter
    k <-    list(rep(0.1, n), rep(0.1, n))
    
    # combine parameters into list 
    para <- list(c = c, opt = opt, k = k)

#------------------#
# 3. Simulation ----
#------------------# 
    
    # This function simultes the abundances. For more details see the collected_functions script 
    out <- sim.re.md(n.grad = 2, re.type = list(rep(3, n), rep(3, n)), parameter = para, grad = env) 

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