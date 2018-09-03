                            ### -------------------------- ###
                            # --- Example GLMmv analysis --- # 
                            ### -------------------------- ###


# This scirpt shows an exemplary GLMmv analysis of the class 1 Uni-Uni model.
    
#-------------#
# -- SETUP ----
#-------------#
    
    # Load packages 
    library(mvabund)
    
    # Set working directory to folder with simulation scripts.
    setwd("~/Uni/Master Thesis/01_R/Geputzt/01_simulation/Communities /")
    
    # Load Model
    source("01_Uni-Uni_Class1.R")

#---------------------#
# -- PREPARE DATA ----
#---------------------#    
    
    # Remove emtpy sites 
    index <- apply(out.samp,1,sum)>0
    
    out.samp <- out.samp[index,]
    envn <-envn[index,]

    # transform abundace data to mvabund-type object
    com_uni <- mvabund(out.samp)
    
    # Plot showing the mean-variance relationship of the data
    meanvar.plot(com_uni, xlab = "Mean", ylab = "Variance")
    
#----------------------#
# -- GLMmv ANALYSIS ----
#----------------------#
        
    # call GLMmvs with different residual distributions
    mod_nb  <- manyglm(com_uni ~ ., data = envn, family = "negative.binomial")
    mod_pos <- manyglm(com_uni ~ ., data = envn, family = "poisson")
    mod_gaus <- manylm(com_uni ~ ., data = envn)

    
#---------------#
# -- RESULTS ----
#---------------#
    
    # Calculate corresponding AICs - one value for each species 
    AIC_table <- data.frame(AIC_NB = AIC(mod_nb),
                                AIC_P = AIC(mod_pos),
                                AIC_G = AIC(mod_gaus))

    # Residual plots 
    plot(mod_nb, which = 1)
    plot(mod_pos, which = 1)
    plot(mod_gaus, which = 1)
    plot(mod_nb, which = 2)
    plot(mod_pos, which = 2)
    plot(mod_gaus, which = 2)
    plot(mod_nb, which = 3)
    plot(mod_pos, which = 3)
    plot(mod_gaus, which = 3)

    
    # Anova analysis 
    an.res <- anova.manyglm(mod_nb, p.uni="adjusted", nBoot = 500)