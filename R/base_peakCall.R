
# calculate posterior probabilities (1 signal model)

.getPH_1S <- function( margDensity, pi0 )
{
    # take parameters
    
    MDZ0 <- margDensity$MDZ0
    MDZ1 <- margDensity$MDZ1
    
    
    # calculate posterior probabilities
    
    denom <- MDZ0*pi0 + MDZ1*(1-pi0)
    pH0 = MDZ0 * pi0 / denom  
    pH1 = MDZ1 * (1-pi0) / denom

    if( length(which(is.na(pH0))) > 0 )
    {
        pH0[ is.na(pH0) ] = 1
        pH1[ is.na(pH0) ] = 0
    }
    post_prob = list( pH0=pH0, pH1=pH1 )
    
    return( post_prob )
}


# calculate posterior probabilities (2 signal model)

.getPH_2S <- function( margDensity, pi0, p1 )
{
    # take parameters
    
    MDZ0 <- margDensity$MDZ0
    MDZ1 <- margDensity$MDZ1
    MDZ2 <- margDensity$MDZ2
    
    
    # calculate posterior probabilities
    
    denom <- MDZ0*pi0 + (1-pi0) * ( MDZ1*p1 + MDZ2*(1-p1) )
    pH0 = MDZ0 * pi0 / denom  
    pH1 = MDZ1 * (1-pi0) * p1 / denom
    pH2 = MDZ2 * (1-pi0) * (1-p1) / denom
    

    if( length(which(is.na(pH0))) > 0 )
    {
        pH0[ is.na(pH0) ] = 1
        pH1[ is.na(pH0) ] = 0
        pH2[ is.na(pH0) ] = 0
    }
    post_prob = list( pH0=pH0, pH1=pH1, pH2=pH2 )
    
    return( post_prob )
}


# claim peaks

.peakCall <- function( postProb, dataSet, FDR, maxgap=200, minsize=50, thres=10, analysisType )
{   
    #library(IRanges)
    
    # FDR = direct posterior probability approach (Newton et al., 2004, Biostatistics)
    
    # determine peaks

    betapH <- postProb$pH0
    betapH_s <- sort(betapH)
    sbetapH <- cumsum(betapH_s) / c(1:length(betapH))		# expected rate of false discoveries
    id <- which( sbetapH <= FDR )
    
    if(length(id)>0)
    {
        #################### if peaks exist
        
        coord <- dataSet$coord
        Y <- dataSet$Y
        switch( analysisType,
            OS = {
                M <- dataSet$M
                GC <- dataSet$GC
            },
            TS = {
                X <- dataSet$X
                M <- dataSet$M
                GC <- dataSet$GC
            },
            IO = {
                X <- dataSet$X
            }
        ) 
                
        
        # threshold peaks by min tag count & determine peaks
        
        cutoff <- betapH_s[max(id)]    
        bd_bin <- rep( 0, length(Y) )
        bd_bin[ betapH<=cutoff & Y>thres ] <- 1
        
        if ( length(which( betapH<=cutoff & Y>thres )) > 0 ) {
            # if we still have peaks after tag count thresholding
             
            # empirical FDR
            
            empFDR <- sum(betapH[which(betapH<=cutoff)]) / length(which(betapH<=cutoff))
            empFDR_thres <- sum(betapH[which( betapH<=cutoff & Y>thres )]) / 
                length(which( betapH<=cutoff & Y>thres ))
            #cat( "Info: empirical FDR (before thresholding) = ", 
            #    round(1000*empFDR)/1000, "\n", sep="" )        
            #cat( "Info: empirical FDR (after thresholding) = ", 
            #    round(1000*empFDR_thres)/1000, "\n", sep="" )
            
            
            # generate initial peak list
      
            bd_ID <- which(bd_bin==1)           # initial peak (bin-level)
            indRanges <- IRanges( start=bd_ID, end=bd_ID+1 )
            reducedRanges <- reduce(indRanges)  # merge nearby peaks
            
            binsize <- coord[2] - coord[1]
            peak_start <- coord[ start(reducedRanges) ]
            peak_stop <- coord[ end(reducedRanges)-1 ] + binsize - 1
            
            
            # merge close peaks if distance<=maxgap
            
            coord_org <- IRanges( start=peak_start, end=peak_stop+maxgap )
            coord_merged <- reduce(coord_org)
            peak_start <- start(coord_merged)
            peak_stop <- end(coord_merged) - maxgap
            
            
            # filter peaks smaller than minsize & order by coordinates
            
            peaksize <- peak_stop - peak_start + 1
            filterID = which( peaksize <= minsize )
            if ( length(filterID) > 0 )
            {
                peak_start <- peak_start[ -filterID ]
                peak_stop <- peak_stop[ -filterID ]
            }
            peak_start <- peak_start[ order(peak_start) ]        
            peak_stop <- peak_stop[ order(peak_start) ]
            peaksize = peak_stop - peak_start
            
            
            # calculate additional info
            
            switch( analysisType,
                OS = {        
                    # extract info
                    
                    peak_info <- apply( cbind(peak_start,peak_stop), 1, 
                        function(x) {
                            start_i <- x[1]
                            end_i <- x[2]
                            ind_i <- which( coord>=start_i & coord<=(end_i-binsize+1) )
                            betapH_i <- betapH[ind_i]
                            aveP <- mean(betapH_i)
                            minP <- min(betapH_i)
                            aveChipCount <- mean(Y[ind_i])
                            maxChipCount <- max(Y[ind_i])
                            aveM <- mean(M[ind_i])
                            aveGC <- mean(GC[ind_i]) 
                            return( c( aveP, minP, aveChipCount, maxChipCount, aveM, aveGC ) )
                        }
                    )
                    
                    # combine all
                            
                    final_peakset <- data.frame( peak_start, peak_stop, peaksize, t(peak_info) )
                    colnames(final_peakset) =
                        c('peakStart','peakStop','peakSize','aveP','minP',
                        'aveChipCount','maxChipCount','map','GC')
                },
                TS = {        
                    # extract info
                            
                    peak_info <- apply( cbind(peak_start,peak_stop), 1, 
                        function(x) {
                            start_i <- x[1]
                            end_i <- x[2]
                            ind_i <- which( coord>=start_i & coord<=(end_i-binsize+1) )
                            betapH_i <- betapH[ind_i]
                            aveP <- mean(betapH_i)
                            minP <- min(betapH_i)
                            aveChipCount <- mean(Y[ind_i])
                            maxChipCount <- max(Y[ind_i])
                            aveInputCount <- mean(X[ind_i])
                            nRatio <- sum(Y) / sum(X)
                            aveInputCountScaled <- aveInputCount * nRatio
                            aveLog2Ratio <- mean( log2( (Y[ind_i]+1) / (X[ind_i]*nRatio+1) ) )
                            aveM <- mean(M[ind_i])
                            aveGC <- mean(GC[ind_i]) 
                            return( c( aveP, minP, aveChipCount, maxChipCount, 
                                aveInputCount, aveInputCountScaled, aveLog2Ratio, aveM, aveGC ) )
                        }
                    )
                    
                    # combine all
                            
                    final_peakset <- data.frame( peak_start, peak_stop, peaksize, t(peak_info) )
                    colnames(final_peakset) =
                    c('peakStart','peakStop','peakSize','aveP','minP',
                    'aveChipCount','maxChipCount',
                    'aveInputCount','aveInputCountScaled','aveLog2Ratio','map','GC')
                },
                IO = {        
                    # extract info
                            
                    peak_info <- apply( cbind(peak_start,peak_stop), 1, 
                        function(x) {
                            start_i <- x[1]
                            end_i <- x[2]
                            ind_i <- which( coord>=start_i & coord<=(end_i-binsize+1) )
                            betapH_i <- betapH[ind_i]
                            aveP <- mean(betapH_i)
                            minP <- min(betapH_i)
                            aveChipCount <- mean(Y[ind_i])
                            maxChipCount <- max(Y[ind_i])
                            aveInputCount <- mean(X[ind_i])
                            nRatio <- sum(Y) / sum(X)
                            aveInputCountScaled <- aveInputCount * nRatio
                            aveLog2Ratio <- mean( log2( (Y[ind_i]+1) / (X[ind_i]*nRatio+1) ) )
                            return( c( aveP, minP, aveChipCount, maxChipCount, 
                                aveInputCount, aveInputCountScaled, aveLog2Ratio ) )
                        }
                    )
                    
                    # combine all
                            
                    final_peakset <- data.frame( peak_start, peak_stop, peaksize, t(peak_info) )
                    colnames(final_peakset) =
                    c('peakStart','peakStop','peakSize','aveP','minP',
                    'aveChipCount','maxChipCount',
                    'aveInputCount','aveInputCountScaled','aveLog2Ratio')
                }
            )
            
            
            #################### if peaks exist
            
            return( list( peakSet=final_peakset, bdBin=bd_bin, empFDR=empFDR ) )
        } else {
            # if we lose all peaks after tag count thresholding
             
            # empirical FDR
            
            empFDR <- sum(betapH[which(betapH<=cutoff)]) / length(which(betapH<=cutoff))
            empFDR_thres <- 0
            message( "Info: no peaks remain after thresholding")
            #message( "Info: empirical FDR (before thresholding) = ", round(1000*empFDR)/1000 )        
            #message( "Info: empirical FDR (after thresholding) = ", round(1000*empFDR_thres)/1000 )
            
            Y <- dataSet$Y
            
            return( list( peakSet=NULL, bdBin=rep(0,length(Y)), empFDR=0 ) )
        
        }
    } else
    {
        #################### if there is no peak
        
        Y <- dataSet$Y
        
        return( list( peakSet=NULL, bdBin=rep(0,length(Y)), empFDR=0 ) )  
    }
}
