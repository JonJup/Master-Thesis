                        ### ------------------------- ###
                        # --- Example CCA analysis --- # 
                        ### ------------------------- ###

# This scirpt shows an exemplary CCA analysis of the class 1 Uni-Uni model.

#-------------#
# -- SETUP ----
#-------------#

    # Load packages 
    library(vegan)
    
    
    # Set working directory to folder with simulation scripts.
    setwd("~/Uni/Master Thesis/01_R/Geputzt/01_simulation/Communities /")
    
    # Load Model 
    source("01_Uni-Uni_Class1.R")
    
#---------------------#
# -- PREPARE DATA ----
#---------------------#

    # Remove emtpy sites 
    index    <- apply(out.samp,1,sum)>0
    
    out.samp <- out.samp[index,]
    envn     <- envn    [index,]

    
#--------------------#
# -- CCA ANALYSIS ----
#--------------------#

    mod <- cca(out.samp ~. , envn)
    
    # Run ANOVAs
    anv1 <- anova(mod, by = "terms", permutations=1000, parallel = 4)
    anv2 <- anova(mod, by = "axis" , permutations=1000, parallel = 4)
    
#--------------#
# -- RESULTS ---
#--------------#    
    
    # Model Summary
    mod.sum <- summary(mod)
    
    # Plots 
    par(mfrow = c(1,1))
    plot(mod, scaling = 1, display = c("sp", "lc", "cn", "wa"), type = "n", main = "LC & WA scores - scaling 1"); points(mod, display = "wa", col = "green", cex = 2, pch = 17);points(mod, display = "lc", col = "grey", pch = 16, cex = .6);  text(mod, display = "species", col = "red", cex = 1.5); text(mod,"cn", col = "blue")
    plot(mod, scaling = 2, display = c("sp", "lc", "cn", "wa"), type = "n", main = "LC & WA scores - scaling 2"); points(mod, display = "lc", col = "black", pch = 16, cex = .3); text(mod, display = "species", col = "red", cex = 1.5); text(mod,"cn", col = "blue")

