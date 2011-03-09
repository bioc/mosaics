
.getParamZ0 <- function(Y) {      
    
    # encode ty to return numeric vector
    # Q -> 1, MM -> 2, NA -> 3
    
    u0 <- length(which(Y==0))
    u1 <- length(which(Y==1))
    u2 <- length(which(Y==2))
    mean0 <- mean(Y)
    var0 <- var(Y)
    
    if( u0>0 & u1>0 & u2>0 ) {
        r1 <- u1/u0
        r2 <- u2/u1
        a <- r1/(2*r2-r1)
        b <- 1/(2*r2-r1)-1
        #ty <- 'Q'
        ty <- 1
    
        if(a<=0||b<=0||a==Inf||b==Inf){
            if(var0 > mean0){
                a <- mean0^2/(var0-mean0)
                b <- mean0/(var0-mean0)
                #ty <- 'MM'
                ty <- 2
            }
            if(var0 <= mean0){
                a <- NA
                b <- NA
                #ty <- NA
                ty <- 3
            }
        }
    }else{
        a <- NA
        b <- NA
        #ty <- NA
        ty <- 3
    }
    return( c( a, b, mean0, var0, u0, u1, u2, length(Y), ty ) )
}
