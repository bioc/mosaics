
.getPi0 <- function( mu, a, Y_freq )
{
    # generate simulated tag counts
    
    b <- a/mu
    N <- length(b)
    Y_sim <- rnbinom(N,a,b/(b+1))
    
    
    # calculate pi0
    
    py0z0 <- dnbinom(0,a,b/(b+1)) 
    py1z0 <- dnbinom(1,a,b/(b+1))
    py2z0 <- dnbinom(2,a,b/(b+1)) 
    
    denom_pi0 <- sum(py0z0) + sum(py1z0) + sum(py2z0) 
    num_pi0 <- sum( Y_freq[1:3] )
    pi0 <- min(num_pi0/denom_pi0,1)
    if(pi0==1)
    {
        pi0 <- min(num_pi0/length(which(Y_sim<=2)),0.99)
    }
    
    return( pi0 )
}
