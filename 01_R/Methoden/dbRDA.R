                        ### -------------------------- ###
                        # --- Example dbRDA analysis --- # 
                        ### -------------------------- ###
# This scirpt shows an exemplary dbRDA analysis of the class 1 Uni-Uni model. 
#--------------#
# -- SETUP ----
#--------------#

    # load package 
    library(vegan)
    
    # Set working directory to folder with simulation scripts.
    setwd("~/Uni/Master Thesis/01_R/Geputzt/01_simulation/Communities /")
    
    # Load model 
    source("01_Uni-Uni_Class1.R")
    
    # Remove emtpy sites 
    index <- apply(out.samp,1,sum)>0
    
    out.samp <- out.samp[index,]
    envn <-envn[index,]
    
#----------------------#    
# -- dbRDA ANALYSIS ----  
#----------------------#
    
    # Run dbRDA
    mod <- capscale(out.samp ~ ., data = envn, add = "lingoes", distance = "bray")
    
#--------------#
# -- RESULTS ---
#--------------# 
        
    # Run ANOVA 
    anv1 <- anova.cca(mod, step  = 1000, by = "term", parallel = 4)
    anv2 <- anova.cca(mod, step  = 1000, by = "axis", parallel = 4)
    
    # Triplots 
    plot(mod, scaling = 1, display = c("sp", "lc", "cn", "wa"), type = "n", main = "LC scores - scaling 1"); points(mod, display = "wa", col = "green", cex = 2, pch = 17);points(mod, display = "lc", col = "grey", pch = 16, cex = .6);  text(mod, display = "species", col = "red", cex = 1.5); text(mod,"cn", col = "blue")
    plot(mod, scaling = 2, display = c("sp", "lc", "cn", "wa"), type = "n", main = "LC scores - scaling 2"); points(mod, display = "lc", col = "black", pch = 16, cex = .3); text(mod, display = "species", col = "red", cex = 1.5); text(mod,"cn", col = "blue")



