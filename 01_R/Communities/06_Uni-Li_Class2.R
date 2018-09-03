### --- Uni-Li Class 2 --- ###  

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
    n <- 5
    
    # number of sampling sites per gradient 
    sites <- 10
    
    # number of noise variables 
    nerr <- 5
    
    # env is the grid of x and y gradient along which species abundances will be simulated.     
    env <- expand.grid(1:100,1:100)

    # Now the species' response pararmeters are set: 
    # maximum abundance 
    c = list("G" = list(rep(100,n), # 1. Grad
                        rep(100,n))) # 2. Grad
    
    # Tolerance 
    tol = list("G" = list(rep(7.5,n),
                          rep(7.5,n)))
    
    #Optimum
    opt = list("G" = list(seq(from = 0, to = 100, length.out = 5)))
    
    # Linear coefficient
    lc = list(rep(0.1,n),c(1,1.2,1.4,1.6,1.8))
    
    # combine all parameters into one list 
    para <- list(c=c, tol=tol,opt=opt, lc = lc)

#------------------#
# 3. Simulation ----
#------------------#   
    
    # This function simultes the abundances. For more details see the collected_functions script.
    out <- sim.re.md(n.grad = 2, re.type = list(rep(2,n),rep(1,n)), parameter = para, grad = env)
    
    # Assign names sp1, sp2 ... to species 
    names(out) <- paste("sp", 1:ncol(out), sep = "")

    sample <- take.sample(out, n.x = sites, n.y = sites, grid = env)
    
    # Extract new abundance data  
    out.samp <- sample[["Abundance"]]
    
    # Extract new sample location
    out.grid <- sample[["Sample grid"]]

    # Add noise variables 
    envn <- make.noise(n.grad = nerr, n.cor = .0, grid = out.grid, h = 100)