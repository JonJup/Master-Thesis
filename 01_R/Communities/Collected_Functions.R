### --- COLLECTED_FUNCTIONS --- ### 

# This script holds the function that are needed for the simulations.  

#--------------#
# Complement()
#--------------#

# This function takes a set of variables y and adjusts a variable x so that x
# correlates to y as determined by rho.

complement <- function(y, rho, x) {
    if(!is.matrix(y)) y <- matrix(y, ncol=1)
    if (missing(x)) x <- rnorm(n)
    d <- ncol(y) 
    n <- nrow(y)
    y <- scale(y)
    e <- residuals(lm(x ~ y))
    y.dual <- with(svd(y), (n-1)*u %*% diag(ifelse(d > 0, 1/d, 0)) %*% t(v))
    sigma2 <- c((1 - rho %*% cov(y.dual) %*% rho) / var(e))
    if (sigma2 >= 0) {
        sigma <- sqrt(sigma2) 
        z <- y.dual %*% rho + sigma*e
    } else {
        warning("Correlations are impossible.")
        z <- rep(0, n)
    }
    return(z)
}

#--------------#
# make.noise()
#--------------#

# This functions creates n.grad noise variables. The length of these vectors is determined by the grid argument. 
# The grid that was used for the simulation of abundance data should be supplied here. 
# n.cor sets the correlation between the noise variables. 
# If the variables are supposed to be scaled, the maximum value can be set by h.  

make.noise <- function(n.grad, grid, n.cor, h){
    
    set.seed(1); RHO <- list();  out <- as.matrix(grid)
    for (i in 1:n.grad) RHO[[i]] <- rep(n.cor, i+1)
    for (i in 1:n.grad){
        rv <- rnorm(n = nrow(out), mean = 0, sd = 1)
        rv <- complement(out, rho = RHO[[i]], rv)
        out <- cbind(out,rv)
    }
    for (i in 1:ncol(out)){
        out[,i] <- out[,i] * (h/max(out[,i]))
    }
    out <- as.data.frame(out); env.names  <- paste("env", 1:ncol(grid),sep = "");rand.names <- paste("rand",1:n.grad, sep = "");    names(out) <- c(env.names,rand.names)
    return(out)
}

#----------------#
# sim.re.md()----
#----------------#

# This function simulates the abundaces.
# n.grad: Integer, number of environmental variables, i.e. columns in grad 
# re.type: List, response type of each species on each gradient. 
#   0 = No response; 1 = Unimodal; 2. Linear; 3. Logisitc; 4. Bimodal
# grad: Matrix, environmental variables
# parametr: List, parameters that are required for the specified response types.

sim.re.md <- function(n.grad, re.type, grad, parameter) {
    
    # prepare 
    n.species <- length(re.type[[1]])   # number of species in model 
    resp.mat <-  matrix(ncol = n.species, nrow = nrow(grad)) #create empty matrix to hold responses  
    resp.list <- list(resp.mat,resp.mat) # create list that holds as many instances of resp.mat as their are dimensions (predictors) in the system. Currently only a maximum of two is supported.
    opt.used <-  list("Grad1" = c(),
                      "Grad2" = c())
    
    # unpack parameter-list so the model becomes easier to read. the suffix u stands for unpacked and is added to avoid double naming. 
    # linear part 
    lc.u <- parameter[["lc"]]
    # gaussian
    c.g.u   <- parameter[["c"]][["G"]]
    opt.g.u <- parameter[["opt"]][["G"]]
    tol.g.u <- parameter[["tol"]][["G"]]
    # logistic 
    c.l.u   <- parameter[["c"]][["L"]]
    k.l     <- parameter[["k"]]
    opt.l.u <- parameter[["opt"]][["L"]]
    # binomial 
    c.b.u <- parameter[["c"]][["B"]]
    opt.b.u <- parameter[["opt"]][["B"]]
    tol.b.u <- parameter[["tol"]][["B"]]
    
    # Indices for responses. If the first linear response was calculed n.lin is set to 2. 
    n.lin <- n.gaus <- n.log <- n.bino <- c(1,1)
    
    # start gradient loop. The abundances of each species is calculated for the first gradient, then the second and so on. 
    for(j in 1:n.grad){
        
        # Start species loop. First the abundance of species 1 is calulated for all sites then for the second and so on. 
        for (i in 1:n.species){
            
            ### --- linear response --- ###     
            if (re.type[[j]][i] == 1){
                
                resp.list[[j]][,i] <- grad[,j] * lc.u[[j]][n.lin[j]]
                
                opt.used[j][[1]][i] <- max(grad[,j])
                n.lin[j] <- n.lin[j] + 1
                
                
            }
            
            ### --- gaussian response --- ###
            if (re.type[[j]][i] == 2){
                
                resp.list[[j]][,i] <- 
                    c.g.u[[j]][n.gaus[j]] * exp(-(grad[,j] - opt.g.u[[j]][n.gaus[j]])^2 / 
                                                    (2*tol.g.u[[j]][n.gaus[j]]^2)) 
                
                opt.used[j][[1]][i] <- opt.g.u[[j]][n.gaus[j]]
                n.gaus[j] <- n.gaus[j] + 1
                
                
                
            } 
            
            ### --- logisitc response --- ###
            #parameter 1= max abu, 2 = steigung, 3 = midpoint
            
            if(re.type[[j]][i] == 3){
                
                resp.list[[j]][,i] <- 
                    c.l.u[[j]][n.log[j]] / (1 + exp (-k.l[[j]][n.log[j]] * (grad[,j] - opt.l.u[[j]][n.log[j]]))) 
                
                
                opt.used[j][[1]][i] <- opt.l.u[[j]][n.log[j]]
                n.log[j] <- n.log[j] + 1
            }  
            
            ### --- binomial response --- ### 
            if (re.type[[j]][i] == 4){
                
                resp.list[[j]][,i] <- 
                    c.b.u[[j]][[n.bino[j]]][1] * exp(-(grad[,j] - opt.b.u[[j]][[n.bino[j]]][1])^2 / (2 * tol.b.u[[j]][[n.bino[j]]][1]^2)) +  
                    c.b.u[[j]][[n.bino[j]]][2] * exp(-(grad[,j] - opt.b.u[[j]][[n.bino[j]]][2])^2 / (2 * tol.b.u[[j]][[n.bino[j]]][2]^2))
                
                
                
                opt.used[j][[1]][i] <- paste(opt.b.u[[j]][[n.bino[j]]][1], 
                                             opt.b.u[[j]][[n.bino[j]]][2], sep = "/")
                
                n.bino[j] <- n.bino[j] + 1
            }    
            ### --- no response --- ### 
            if (re.type[[j]][i] == 0){
                
                resp.list[[j]][,i] <- rep(1,nrow(grad))
                opt.used[j][[1]] <- append(opt.used[j][[1]], 0)
            }
        }    
        
        
    }
    
    mult.resp.mat <- resp.list[[1]] * resp.list[[2]]
    mult.resp.mat <- as.data.frame(round(mult.resp.mat))
    return(mult.resp.mat)
}

#----------------#
# take.sample()
#----------------#

# This function samples form the output of sim.re.md(). 
# data: The data that will be sampled 
# n.x: number of new samples along x axis
# n.y: number of new samples along y axis
# grid: grid that was used in sim.re.md() call

take.sample <- function(data, n.x, n.y, grid){
    
    output <- matrix(nrow = n.x * n.y , ncol = ncol(data))
    
    #regular
    sampling.locs.x <- round(seq(from = min(grid[,1]), to = max(grid[,1]), length.out = n.x))
    sampling.locs.y <- round(seq(from = min(grid[,2]), to = max(grid[,2]), length.out = n.y))
    
    sample.grid <- expand.grid(sampling.locs.x, sampling.locs.y)
    
    index <- sample.grid[,1] + (sample.grid[,2] - 1) * length(unique(grid[,1])) 
    
    output <- data[index,]

    sample.grid <- as.data.frame(sample.grid)
    return(list("Abundance" = output, "Sample grid" = sample.grid))
}
