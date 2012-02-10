#########################################################
# Compute marginal density for Y=N+S1+S2
#########################################################

.margDist_1S <- function( mosaicsEst, tagCount, pNfit, k=3 )
{            
    # extract parameters
    
    a <- mosaicsEst@a
    mu_est <- mosaicsEst@muEst
    b_est <- a / mu_est
        
    b <- mosaicsEst@b
    c <- mosaicsEst@c
    
    Yori <- tagCount
    
    
    # use only Y >= k
    
    Y <- Yori - k 
    Y[which(Y<0)] = -1  
    
    
    # round mu

    mu_round <- round(mu_est,2)
    if( length(which(Y<0))>0 ) mu_round[which(Y<0)] <- 0
    mu_round_U <- unique(mu_round)
    #n_mu_U <- length(mu_round_U)
    
    
    # prob of N using rounding mu for prob of S
     
    Ymax <- max(Y)
    #pN = matrix( 0, nrow=n_mu_U, ncol=length(c(0:Ymax)) ) 
    #for( i in 1:n_mu_U )
    #{
    #    b_round_i = a / mu_round_U[i]
    #    pN[i,] = dnbinom( 0:Ymax, a, b_round_i/(b_round_i+1) )
    #}
    
    pN <- pNfit$pN
    
    
    # prob of N (Yori & b_est have same length)
    
    MDZ0 = dnbinom( Yori, size=a, b_est/(b_est+1) )
    
    
    # prob of S (positive only when Y>=k)
    
    pS <- dnbinom( 0:Ymax, b, c/(c+1) )
   
    id_Y_ge_k <- which( Y>=0 )
    #n_Y_ge_k <- length( id_Y_ge_k )
    MDZ1 <- rep( 0, length(Y) )

    #print( 'Now analyzing with foreground with p1(N+S1)+(1-p1)(N+S2)' )
    #print( paste('TOTAL ITERATION=',n.y.pos) )
    
    #for( i in 1:n_Y_ge_k )
    #{                        
    #    id_i <- id_Y_ge_k[i]
    #    Y_i <- Y[ id_i ]
    #    mu_i <- mu_round[ id_i ]     
    #    pN_i <- pN[ mu_round_U==mu_i, 1:(Y_i+1) ]
    #    
    #    MDZ1[ id_i ] <- sum( rev(pS[1:(Y_i+1)]) * pN_i )        
    #    #if( (n.y.pos-i)%%10000 == 0 ) print( paste(n.y.pos-i,'to go') )
    #}
    
    MDZ1[ id_Y_ge_k ] <- conv_1S( y=Y[id_Y_ge_k], mu_round=mu_round[id_Y_ge_k],
    	mu_round_U=mu_round_U, pN=pN, pS=pS )
    	    
    
    # return object

    MD=data.frame( cbind( MDZ0, MDZ1 ) )
    colnames(MD)=c('MDZ0','MDZ1')
    return(MD)   
}
