### --- Uni-Bi Class 3 --- ###  

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
    sites <- 15
    
    # number of noise variables 
    nerr <- 5
    
    # env is the grid of x and y gradient along which species abundances will be simulated.     
    env <- expand.grid(1:100,1:100)

    # Now the species' response pararmeters are set: 
    # Maximal abundances
    c = list(   "G" = list(rep(100,n),rep(100,n)), 
                "B" = list(list(c(100,100),c(100,100), c(100,100), c(100,100), c(100,100),c(100,100),c(100,100), c(100,100), c(100,100)),
                       list(c(100,100),c(100,100), c(100,100), c(100,100), c(100,100),c(100,100),c(100,100), c(100,100), c(100,100)))
        ) 
    # Tolerance 
    tol = list( "G" = list(rep(5,n), rep(5,n)),
                "B" = list(list(c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5)),
                      list(c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5),c(5,5))))
    
    # Optima 
    opt = list("G" = list(rep(c(20,50,80), times = 3), rep(c(20,50,80), times = 3)),
           "B" = list(list(c(10,30), c(10,30), c(10,30), c(40,60), c(40,60), c(40,60),  c(70,90),c(70,90), c(70,90)),
                      list(c(10,30), c(10,30), c(10,30), c(40,60), c(40,60), c(40,60),  c(70,90),c(70,90), c(70,90))))
               
    # combine all parameters into one list 
    para <- list(c=c, tol=tol,opt=opt)

#------------------#
# 3. Simulation ----
#------------------# 
    
    # This function simultes the abundances. For more details see the collected_functions script   
    out <- sim.re.md(n.grad = 2, re.type = list(rep(2,n),rep(4,n)), parameter = para, grad = env)

    names(out) <- paste("sp", 1:ncol(out), sep = "")

    # Sample the big 100x100 dataset
    sample <- take.sample(out, n.x = sites, n.y = sites, grid = env)
    
    # Extract new abundance data
    out.samp <- sample[["Abundance"]]
    
    # Extract new sample location 
    out.grid <- sample[["Sample grid"]]
    
    # Add noise variables 
    envn <- make.noise(n.grad = nerr, n.cor = .0, grid = out.grid, h = 100)