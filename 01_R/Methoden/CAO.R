            ### ------------------------ ###
            # --- Example CAO analysis --- # 
            ### ------------------------ ###

# This scirpt shows an exemplary CAO analysis of the class 1 Uni-Uni model
#-------------#
# -- SETUP ----
#-------------#

    # Load packages 
    library(VGAM)
    library(stringr)
    library(dplyr)

    # Set working directory to folder with simulation scripts.
    setwd("~/Uni/Master Thesis/01_R/Geputzt/01_simulation/Communities /")
    # Load Model 
    source("01_Uni-Uni_Class1.R")
    
#--------------------#    
# -- PREPARE DATA ----    
#--------------------#    
    
    # Remove emtpy sites 
    index <- apply(out.samp,1,sum)>0
    
    out.samp <- out.samp[index,]
    envn <-envn[index,]

    # combined data frame 
    out.comb <- cbind(out.samp, envn)

#--------------------#
# -- CAO ANALYSIS ----
#--------------------#    
    
    # run model
    mod <- cao(cbind(sp1, sp2, sp3, sp4, sp5,sp6, sp7, sp8, sp9) ~ env1 + env2 + rand1,
               data = out.comb, family =  poissonff, Rank = 1, df1.nl = 1.5, Bestof = 10, Crow1positive = T) 

#---------------#
# -- RESULTS ----
#---------------# 
       
    # summary of results 
    summary.rrvgam(mod)
    
    # history i.e. deviance of seperate model runs. 
    sort(deviance(mod, history = TRUE))
   
    # Plots
     par(mfrow=c(3,3))
    plot(mod, lcol = "blue", lwd = 2, ylim = c(-5,5)) 
   
    par(mfrow=c(1,1))
    S = ncol(mod@y) # Number of species 
    clr = (1:(S + 1))[-7] # Omit yellow color
    persp(mod, col = clr, lwd = 2, label = TRUE)
    
