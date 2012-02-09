
# calculate P(Y|Z=1) for the current parameters

.calcPN <- function( Y, k, a, mu_est ) 
{             
    b_est <- a / mu_est
    
    # parameters for Y>=k
    
    id_geqk <- which(Y>=k)
    
    Yori <- Y[id_geqk]
    b_est <- b_est[id_geqk]
    mu_est <- mu_est[id_geqk]
    
    # P( Y | Z0 )
    
    PYZ0 <- dnbinom( Yori, a, b_est / (b_est+1) )
    
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
    
    return( list( k=k, pN=pN, mu_round=mu_round, mu_round_U=mu_round_U, PYZ0=PYZ0 ) )
}
