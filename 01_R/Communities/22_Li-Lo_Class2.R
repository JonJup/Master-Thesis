### --- Li-Lo Class 2 --- ###  

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
    # linear parameter 
    lc = list(seq(from = 0.1, to = 0.1, length.out = n ),
              seq(from = 0.1, to = 0.1, length.out = n ))
    
    # maximal abundances 
    c <-list("L" = list(rep(100,n),rep(100,n)))
    
    # optima 
    opt <-  list("L" = list( 
        c(15,30,45,60,75),
        c(15,30,45,60,75)))
    
    # logistic parameter
    k <-    list(
        rep(0.05,n),
        rep(0.05,n)
    )

    # combine parameters into a list
    para <- list(lc = lc, k = k, c = c, opt = opt )

#------------------#
# 3. Simulation ----
#------------------# 
    
    # This function simultes the abundances. For more details see the collected_functions script     
    out <- sim.re.md(n.grad = 2, re.type = list(rep(1,n),rep(3,n)), parameter = para, grad = env)
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