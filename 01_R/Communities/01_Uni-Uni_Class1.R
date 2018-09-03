### --- Uni-Uni Class 1 --- ###  

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
# maximum abundance 
c = list("G" = list(rep(100,n),  # 1. Grad
                    rep(100,n))) # 2. Grad
# tolerance 
tol = list("G" = list(rep(7.5,n),
                      rep(7.5,n)))
# optima 
opt = list("G" = list(rep(c(20,50,80), times = 3),
                      rep(c(20,50,80), each = 3)))

# combine all parameters into one list 
para <- list(c=c, tol=tol,opt=opt)


#------------------#
# 3. Simulation ----
#------------------#   

# This function simultes the abundances. For more details see the collected_functions script.
out <- sim.re.md(n.grad = 2, re.type = list(rep(2,n),rep(2,n)), parameter = para, grad = env)
out <- as.data.frame(out)
# Assign names sp1, sp2 ... to species 
colnames(out) <- paste("sp", 1:ncol(out), sep = "")

# Sample the big 100x100 dataset
    sample <- take.sample(out, n.x = sites, n.y = sites, grid = env)
    
    # Extract new abundance data 
    out.samp <- sample[["Abundance"]] 
    
    # Extract new sample location 
    out.grid <- sample[["Sample grid"]]
    
    # Add noise variables 
    envn <- make.noise(n.grad = nerr, n.cor = .0, grid = out.grid, h = 100)
    