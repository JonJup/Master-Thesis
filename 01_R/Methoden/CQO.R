# --- Example CQO analysis --- # 

# ------------#
# -- SETUP ----
# ------------#

    # load packages 
    library(VGAM)
    library(dplyr)
    library(stringr)
    
    # Set working directory to folder with simulation scripts.
    setwd("~/Uni/Master Thesis/01_R/Geputzt/01_simulation/Communities /")
    
    # load model 
    source("01_Uni-Uni_Class1.R")
    
#--------------------#
# -- PREPARE DATA ----
#--------------------#
    
    # Remove emtpy sites 
    index    <- apply(out.samp,1,sum)>0
    
    out.samp <- out.samp[index,]
    envn     <- envn    [index,]

    # number of species and variables 
    nspec <- ncol(out.samp)
    nvar  <- ncol(envn)
    
    # combined data frame 
    out.comb <- cbind(out.samp, envn)

    
#--------------------#
# -- CQA ANALYSIS ----
#--------------------#
   
    
    # Run rank-1 models with all tolerance settings 
    cqo.1.EqTol   <- cqo(cbind(sp1, sp2, sp3, sp4, sp5,sp6, sp7, sp8, sp9) ~ env1 + env2 + rand1,
                           data = out.comb, family =  poissonff,   Rank = 1, df1.nl = 1.5, Bestof = 10, 
                           Crow1positive = T, eq.tolerance = T ) 
    cqo.1.UneqTol <- cqo(cbind(sp1, sp2, sp3, sp4, sp5,sp6, sp7, sp8, sp9) ~ env1 + env2 + rand1,
                           data = out.comb, family =  poissonff,   Rank = 1, df1.nl = 1.5, Bestof = 10,
                           Crow1positive = T, eq.tolerance = F ) 
    cqo.1.ITol    <- cqo(cbind(sp1, sp2, sp3, sp4, sp5,sp6, sp7, sp8, sp9) ~ env1 + env2 + rand1,
                           data = out.comb, family =  poissonff,   Rank = 1, df1.nl = 1.5, Bestof = 10,
                           Crow1positive = T, I.tolerances = T ) 
    
    # Run rank-2 models with all tolerance settings 
    cqo.2.EqTol   <- cqo(cbind(sp1, sp2, sp3, sp4, sp5,sp6, sp7, sp8, sp9) ~ env1 + env2 + rand1,
                         data = out.comb, family =  poissonff,   Rank = 2, df1.nl = 1.5, Bestof = 10, 
                         Crow1positive = T, eq.tolerance = T ) 
    cqo.2.UneqTol <- cqo(cbind(sp1, sp2, sp3, sp4, sp5,sp6, sp7, sp8, sp9) ~ env1 + env2 + rand1,
                         data = out.comb, family =  poissonff,   Rank = 2, df1.nl = 1.5, Bestof = 10,
                         Crow1positive = T, eq.tolerance = F ) 
    cqo.2.ITol    <- cqo(cbind(sp1, sp2, sp3, sp4, sp5,sp6, sp7, sp8, sp9) ~ env1 + env2 + rand1,
                         data = out.comb, family =  poissonff,   Rank = 2, df1.nl = 1.5, Bestof = 10,
                         Crow1positive = T, I.tolerances = T ) 
  
 
#--------------#
# -- RESULTS ----
#--------------#      
       
    # Which fits better rank one or two and which tolerance setting fits the best?
    AIC_table <- c(
                            AIC_1_e = AIC(cqo.1.EqTol),
                            AIC_1_u = AIC(cqo.1.UneqTol),
                            AIC_1_i = AIC(cqo.1.ITol),
                            
                            AIC_2_e = AIC(cqo.2.EqTol),
                            AIC_2_u = AIC(cqo.2.UneqTol),
                            AIC_2_i = AIC(cqo.2.ITol)
                        )
    
    sort(AIC_table)
    
    # Summary
    (Summary.2.EqTol   <- summary.qrrvglm(cqo.2.EqTol))
    
    # Histories 
    (History.2.EqTol       <- sort(deviance(cqo.2.EqTol, history = TRUE)))


    # Plots 
    par(mfrow = c(1,1))
    S <- ncol(depvar(cqo.2.EqTol))
    clr <- (1: (S + 1))[-7] # omit yellow 
    persp(cqo.2.EqTol, label = T, col = clr, main = "Poisson - EqTol")
    
    lvplot(cqo.2.EqTol, label = T, ellipse = T, C = T, ccol = "red", sites = TRUE, scol = "gray80",
         pcol = "blue", pch = 16,  chull = TRUE, main = "(a)")
    
