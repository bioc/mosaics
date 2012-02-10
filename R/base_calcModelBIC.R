
.calcModelBIC <- function( fitZ1, Y, pNfit, k=3, model="2S", type="BIC", npar )
{       
    if ( model=="2S" )
    {
        # parameter estimates
        
        mu_est <- fitZ1$muEst
        a <- fitZ1$a
        pi0 <- fitZ1$pi0
        b_est <- a / mu_est
        p1 <- fitZ1$p1
        b1 <- fitZ1$b1
        c1 <- fitZ1$c1
        b2 <- fitZ1$b2
        c2 <- fitZ1$c2
        
        # choose Y>=k
        
        id_geqk <- which(Y>=k)      
        Y_ori <- Y[id_geqk]
        b_est <- b_est[id_geqk]
        mu_est <- mu_est[id_geqk]
        
        # calculate log likelihood
        
        #PYZ0 <- dnbinom( Y_ori, a, b_est/(b_est+1) )    
        PYZ0 <- pNfit$PYZ0 
        #pNfit <- .calcPN( Y_ori, k, a, mu_est ) 
        PYZ1 <- .margDistZ1_2S( Y_ori, pNfit, b1, c1, b2, c2 )
        PYZ1G1 <- PYZ1$MDG1
        PYZ1G2 <- PYZ1$MDG2    
        
        logLik0 <- log( pi0*PYZ0 + (1-pi0)*( p1*PYZ1G1 + (1-p1)*PYZ1G2) )
        logLik1 <- sum( logLik0[!is.na(logLik0)] )
    } else if ( model=="1S" )
    {
        # parameter estimates
        
        mu_est <- fitZ1$muEst
        a <- fitZ1$a
        pi0 <- fitZ1$pi0
        b_est <- a / mu_est
        b <- fitZ1$b
        c <- fitZ1$c
        
        # choose Y>=k
        
        id_geqk <- which(Y>=k)      
        Y_ori <- Y[id_geqk]
        b_est <- b_est[id_geqk]
        mu_est <- mu_est[id_geqk]
        
        # calculate log likelihood
            
        #PYZ0 <- dnbinom( Y_ori, a, b_est/(b_est+1) )
        PYZ0 <- pNfit$PYZ0
        #pNfit <- .calcPN( Y_ori, k, a, mu_est ) 
        PYZ1 <- .margDistZ1_1S( Y_ori, pNfit, b, c )
            
        logLik0 <- log( pi0*PYZ0 + (1-pi0)*PYZ1 )
        logLik1 <- sum( logLik0[!is.na(logLik0)] )    
    }
    
    # calculate BIC
    
    switch( type,
        "AIC" = {
            penalty <- 2            
        },
        "BIC" = {
            n <- length(Y_ori)
            penalty <- log(n)
        }
    )
    val <- -2 * logLik1 + penalty * npar
    
    return(val)
}
