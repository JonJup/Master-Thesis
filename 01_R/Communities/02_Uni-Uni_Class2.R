### --- Uni-Uni Class 2 --- ###  

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
    
    # Assign names sp1, sp2 ... to species 
    names(out) <- paste("sp", 1:ncol(out), sep = "")

    # Sample the big 100x100 dataset
    sample <- take.sample(out, n.x = sites, n.y = sites, grid = env)
    
    # Extract new abundance data 
    out.samp <- sample[["Abundance"]]
    
    # Extract new sample location 
    out.grid <- sample[["Sample grid"]]
    
    # Add noise variables 
    envn <- make.noise(n.grad = nerr, n.cor = .0, grid = out.grid, h = 100)
    
#------------------#
#   4. Plot --------
#------------------# 
 
library(lattice)
names <- paste("sp",1:ncol(out.samp), sep = "")
names.env <- append(names, c("env1","env2"))
wire <- cbind(out.samp,out.grid)
colnames(wire) <- names.env
trellis.par.set("axis.line",list(col=NA,lty=1,lwd=1))
trellis.par.set(name = "fontsize", list(text =14, points = 500))
plot.wire <- wireframe(as.formula(paste(paste(names,collapse="+"), "~ env1 * env2",sep = "")),
                       data=wire, xlab="env1", ylab="env2", zlab=list("",rot = 0), aspect = c(100/100,1),
screen = list(z = 20, x = -65, y = 0), par.box = list(col=NA), colorkey = F
)
plot.wire



