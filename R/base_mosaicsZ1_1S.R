
#########################################################
# Z1 with 1 signal component
#########################################################

.mosaicsZ1_1S <- function( MOSAiCS_Z0, Y, k=3 )
{
    ##################################################################
    ### initialization of the main EM
    ##################################################################
    
    a <- MOSAiCS_Z0$a
    #mu_est <- get_mu_est(M, GC, MOSAiCS_Z0, trunc_GC=0.55)
    mu_est <- MOSAiCS_Z0$muEst
    b_est <- a / mu_est
    pi0 <- MOSAiCS_Z0$pi0
    tab_Y <- MOSAiCS_Z0$Y_freq
    Y_u <- MOSAiCS_Z0$Y_val
    
    
    # Initialization for initial EM iteration
    
    #Eyz0 <- Y_Z1_tmp <- rep(0,length(Y_u))
    Eyz0 <- rep(0,length(Y_u))                  # N_Z0(y)
    Eyz0[1] <- round(sum(dnbinom(Y_u[1],a,b_est/(b_est+1))))

    i <- 2
    while( Eyz0[i-1]>0 & i<=length(Y_u) )
    {
        Eyz0[i] <- round(sum(dnbinom(Y_u[i],a,b_est/(b_est+1))))
        i <- i + 1      
    }

    Y_Z1_tmp <- tab_Y - pi0*Eyz0                # (1-pi0) * N_Z1(y) = N(y) - pi0 * N_Z0(y)
    Y_Z1_tmp[which(Y_Z1_tmp<0)] <- 0

    Y_bd <- cbind( Y_u[-c(1:3)]-k, Y_Z1_tmp[-c(1:3)] ) ### already adjusted for k
    
    if ( sum(Y_bd[,2])==0 ) {
        # nothing left after background are excluded from the data
        # risk level = extremely high
        stop( "over-estimation of background detected. Please tune the parameters!" )
    } else if ( Y_bd[which.max(Y_bd[,2]),2] / sum(Y_bd[,2]) >= 0.90 ) {
    		# almost nothing left after background are excluded from the data
        # risk level = very high
        stop( "over-estimation of background detected. Please tune the parameters!" )
    } else if ( Y_bd[which.max(Y_bd[,2]),2] / sum(Y_bd[,2]) >= 0.80 ) {
    		# still not much left after background are excluded from the data
    		# rick level = medium high
    		warning( "there is some possibility of background over-estimation." )
    		warning( "parameter tuning might provide better fits." )
    } else {
        Y_bd_all <- rep(Y_bd[,1],Y_bd[,2])
    }        
    
    # initial mu & var

    mu1_init <- mean(Y_bd_all)
    var1_init <- var(Y_bd_all)
    
    
    # calculate initial EN and varN, based on initial mu & var

    EN <- mean(mu_est)
    varN <- EN*( 1 + EN/a ) 

    if( mu1_init < EN )
    {
        EN <- mu1_init - 0.1
    }

    if( var1_init < varN )
    {
        varN <- var1_init - 0.15
    }
    
    
    # initialize b & c

    if( var1_init - varN - mu1_init + EN <= 0 )
    {
        b_init <- (mu1_init - EN)^2 / 0.1
        c_init <- (mu1_init - EN) / 0.1   
    } else
    {
        b_init <- (mu1_init - EN)^2 / (var1_init - varN - mu1_init + EN)
        c_init <- (mu1_init - EN) / (var1_init - varN - mu1_init + EN)    
    }  


    ##################################################################
    ### The main EM calculation: iterate till convergence
    ##################################################################
    
    # calculate EN and vanN
    # (redefine EN after heuristic EN adjustment for initialization)

    EN <- mean(mu_est)
    varN <- EN*(1 + EN/a) 
    
    
    # parameters for Y>=k
    
    id_geqk <- which(Y>=k)
    Y_ori <- Y[id_geqk]
    Y_Z1 <- Y - k               # adjusted Y
    Y_Z1 <- Y_Z1[id_geqk]
    #M_Z1 <- M[id_geqk]
    #GC_Z1 <- GC[id_geqk]
    b_est_Z1 <- b_est[id_geqk]
    mu_est_Z1 <- mu_est[id_geqk]
    
    
    # Initialization of the main EM calculation
    
    PYZ0 <- dnbinom( Y_ori, a, b_est_Z1/(b_est_Z1+1) )
    pNfit <- .calcPN( Y_ori, k, a, mu_est_Z1 ) 
    PYZ1 <- .margDistZ1_1S( Y_ori, pNfit, b_init, c_init )

    b_iter <- b_init
    c_iter <- c_init
    
    
    # main EM iteration

    eps <- 1e-6
    #logLik <- -Inf
    #logLik <- c(logLik, sum(log(pi0*PYZ0 + (1-pi0)*PYZ1)))
    logLik1 <- sum( log( pi0*PYZ0 + (1-pi0)*PYZ1 ) )
    logLik <- c( -Inf, logLik1 )
    iter <- 2
    
    while( abs(logLik[iter]-logLik[iter-1])>eps & iter < 10 )
    {
        # E-step
        
        Z_latent <- (1 - pi0)*PYZ1 / ( pi0*PYZ0 + (1 - pi0)*PYZ1 )
        
        # M-step
        
        mu1 <- sum(Z_latent*Y_Z1) / sum(Z_latent)
        var1 <- sum(Z_latent*(Y_Z1 - mu1)^2) / sum(Z_latent)

        b <- (mu1 - EN)^2 / (var1 - varN - mu1 + EN)
        c <- (mu1 - EN) / (var1 - varN - mu1 + EN)
        
        # stop iteration if assumptions are not satisfied
        
        if ( b<0 | c<0 )
        {
            b <- b_iter[(iter-1)]
            c <- c_iter[(iter-1)]
            break
        }
        
        # calculate P(Y|Z=1)
        
        PYZ1 <- .margDistZ1_1S( Y_ori, pNfit, b, c )
        
        # update iteration
        
        #logLik <- c(logLik, sum(log(pi0*PYZ0 + (1-pi0)*PYZ1)))
        logLik_t <- sum( log( pi0*PYZ0 + (1-pi0)*PYZ1 ) )
        logLik <- c( logLik, logLik_t )
        b_iter <- c(b_iter, b)
        c_iter <- c(c_iter, c)

        iter <- iter + 1
    }
    

    #return(list(M_u = MOSAiCS_Z0$M_u, GC_u = MOSAiCS_Z0$GC_u, a_u = MOSAiCS_Z0$a_u, b_u = MOSAiCS_Z0$b_u, mean0_u = MOSAiCS_Z0$mean0_u, var0_u = MOSAiCS_Z0$var0_u, u0_u = MOSAiCS_Z0$u0_u, u1_u = MOSAiCS_Z0$u1_u, u2_u = MOSAiCS_Z0$u2_u, n_u = MOSAiCS_Z0$n_u, ty_u = MOSAiCS_Z0$ty_u, Y_val = MOSAiCS_Z0$Y_val, Y_freq = MOSAiCS_Z0$Y_freq, pi0 = MOSAiCS_Z0$pi0, a = MOSAiCS_Z0$a, beta0_strata = MOSAiCS_Z0$beta0_strata, betaM_strata = MOSAiCS_Z0$betaM_strata, betaM2_strata = MOSAiCS_Z0$betaM2_strata, betaGC_strata = MOSAiCS_Z0$betaGC_strata, b = b, c = c, b_init = b_init, c_init = c_init, logLik = logLik, b_iter = b_iter, c_iter = c_iter))
    
    #return( list( 
    #    M_u = MOSAiCS_Z0$M_u, GC_u = MOSAiCS_Z0$GC_u, a_u = MOSAiCS_Z0$a_u, b_u = MOSAiCS_Z0$b_u,
    #    mean0_u = MOSAiCS_Z0$mean0_u, var0_u = MOSAiCS_Z0$var0_u, n_u = MOSAiCS_Z0$n_u, ty_u = MOSAiCS_Z0$ty_u,
    #    Y_val = MOSAiCS_Z0$Y_val, Y_freq = MOSAiCS_Z0$Y_freq,
    #    pi0 = MOSAiCS_Z0$pi0, a = MOSAiCS_Z0$a, b = b, c = c,
    #    mu_est = mu_est ) )
    
    return( list( 
        pi0 = MOSAiCS_Z0$pi0, a = MOSAiCS_Z0$a, muEst = MOSAiCS_Z0$muEst,
        b = b, c = c ) )
}


# calculate P(Y|Z=1) for the current parameters

.calcPN <- function( Yori, k, a, mu_est ) 
{    
    # process Y
    
    Y <- Yori - k        # use only Y >= k
    if(length(which(Y<0))>0 ) Y[which(Y<0)] <- -1       
    
    # round mu

    mu_round <- round(mu_est,2)
    if(length(which(Y<0))>0 ) mu_round[which(Y<0)] <- 0
    mu_round_U <- unique(mu_round)  
    
    # prob of N (using rounding mu for prob of S)
     
    Ymax <- max(Y)
    pN <- apply( as.matrix(mu_round_U), 1, 
        function(x) {
            b_round <- a / x
            return( dnbinom( 0:Ymax, a, b_round/(b_round+1) ) )
        }
    )
    pN <- t(pN)
    
    return( list( k=k, pN=pN, mu_round=mu_round, mu_round_U=mu_round_U ) )
}

.margDistZ1_1S <- function( Yori, pNfit, b, c )
{     
    k <- pNfit$k
    pN <- pNfit$pN
    mu_round <- pNfit$mu_round
    mu_round_U <- pNfit$mu_round_U
    
    # process Y
    
    Y <- Yori - k        # use only Y >= k
    if(length(which(Y<0))>0 ) Y[which(Y<0)] <- -1 
    Ymax <- max(Y)
    ind_ge_k <- which(Y>=0)
    
    # prob of S
    
    pS <- dnbinom( 0:Ymax, b, c/(c+1) )   
    MDZ1 <- rep( 0, length(Y) )            
    MDZ1[ ind_ge_k ] <- apply( cbind( Y, mu_round )[ ind_ge_k, ], 1, 
        function(x) {
            Y.i <- x[1]
            mu.i <- x[2]
            pN.i <- pN[ mu_round_U==mu.i, 1:(Y.i+1) ]     
            return( sum( rev(pS[1:(Y.i+1)]) * pN.i ) )
        }
    )
    
    return(MDZ1)   
}
