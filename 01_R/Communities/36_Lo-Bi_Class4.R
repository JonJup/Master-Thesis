### --- Lo-Bi Class 4 --- ###  

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
sites <- 5

# number of noise variables 
nerr <- 5

# env is the grid of x and y gradient along which species abundances will be simulated.     
env <- expand.grid(1:100,1:100)

# Now the species' response pararmeters are set:
# maximal abundance 
c <-    list("L" = list( 
    rep(100,n),
    rep(100,n)
                        ),
                 "B" = list(
                        list(c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100)),
                        list(c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100),c(100,100))
                 )
            )
    
    # Optima                     
    opt <-  list("L" = list( 
                        c(20, 50, 80, 20, 50, 80, 20, 50, 80),
                        c(20, 20, 20, 50, 50, 50, 80, 80 ,80)
                        ),
                "B" = list(
                        list(c(5,25),c(5,25),c(5,25),c(35,55),c(35,55),c(35,55),c(75,95),c(75,95),c(75,95)),
                        list(c(5,25),c(5,25),c(5,25),c(35,55),c(35,55),c(35,55),c(75,95),c(75,95),c(75,95))
                ) 
               )
      
    # logisitc parameter          
    k <-    list(
                        rep(0.1, n),
                        rep(0.1, n)
                        )
    
    # Tolerances 
    tol <- list ( "B" = list ( list(c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6)),
                               list(c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6),c(6,6))))
    
    # combine Parameters into one list 
    para <- list(c = c, opt = opt, k = k, tol = tol)

    #------------------#
    # 3. Simulation ----
    #------------------# 
    
    # This function simultes the abundances. For more details see the collected_functions script 
    out <- sim.re.md(n.grad = 2, re.type = list(rep(3, 9), rep(4, 9)), parameter = para, grad = env) 
    
    out <- as.data.frame(out)
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