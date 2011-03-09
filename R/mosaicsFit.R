
# fit MOSAiCS

setMethod(
    f="mosaicsFit",
    signature="BinData",
    definition=function( object, analysisType=NULL, k=3,
        meanThres=NA, s=2, d=0.25, truncProb=0.9999 )
    {
        # Note: users can tune parameters only regarding MOSAiCS model fitting.
        # Note: tuning adaptive griding parameters is not supported yet.
        
        if ( is.null(analysisType) ) {
        		# if "analysisType" is not specified, "analysisType" is determined by dataset
        		
        		if ( length(object@input)==0 ) {
        				# we don't have input: we should have both of M and GC
        				
        				if ( length(object@mappability)>0 & length(object@gcContent)>0 ) {
            				#cat( "Info: one-sample analysis.\n" )
        						analysisType <- "OS"
        				} else {
            				stop( "any of control data, mappability, or GC content does not exist. Cannot proceed!\n" )
        				} 
        		} else {
        				# we have input: TS if we have both M & GC; IO otherwise. 
        				
        				if ( length(object@mappability)>0 & length(object@gcContent)>0 ) {
            				#cat( "Info: two-sample analysis (with mappability & GC content).\n" )
        						analysisType <- "TS"
        				} else {
            				#cat( "Info: two-sample analysis (Input only).\n" )
        						analysisType <- "IO"
        				}
        		}
        } else {
						# if "analysisType" is specified, check its validity
						        
		        # error treatment: Input-only analysis is impossible if input data does not exist
		        
		        if ( analysisType=="IO" & length(object@input)==0 )
		        {
		            message( "Info: two-sample analysis (Input only)." )
		            stop( "control data does not exist. Cannot proceed!\n" )
		        }
		        
		        # error treatment: If M or GC does not exist, TS analysis is not available
		        
		        if ( analysisType=="OS" & ( length(object@mappability)==0 | length(object@gcContent)==0 ) )
		        {
		            message( "Info: one-sample analysis." )
		            stop( "mappability or GC content does not exist. Cannot proceed!\n" )
		        }
		        
		        # error treatment: If input data does not exist, TS analysis is not available
		        
		        if ( analysisType=="TS" & length(object@input)==0 )
		        {
		            message( "Info: two-sample analysis (with mappability & GC content)." )
		            message( "Info: control data does not exist." )
		            message( "Info: one-sample analysis will be implemented instead." )
		            analysisType <- "OS"
		        }
		        
		        # error treatment: If M or GC does not exist, TS analysis is not available
		        
		        if ( analysisType=="TS" & ( length(object@mappability)==0 | length(object@gcContent)==0 ) )
		        {
		            message( "Info: two-sample analysis (with mappability & GC content)." )
		            message( "Info: mappability or GC content data does not exist." )
		            message( "Info: two-sample analysis (Input only) will be implemented instead." )
		            analysisType <- "IO"
		        }
		    }
        
        # default meanThres for each of "OS" & "TS"
        
        if ( is.na(meanThres) )
        {
            switch( analysisType,
                OS = {
                    meanThres <- 0
                },
                TS = {
                    meanThres <- 1
                },
                IO = {
                    meanThres <- NA     # meanThres is not used for analysisType=="IO"
                }
            )
        }
        
        # MOSAiCS model fit
        
        switch( analysisType,
            OS = {
                # one-sample analysis
                
                message( "Info: one-sample analysis." )
                fit <- .mosaicsFit_OS( object, k=k, meanThres=meanThres )
            },
            TS = {
                # two-sample analysis (with M & GC)
                
                message( "Info: two-sample analysis (with mappability & GC content)." )
                fit <- .mosaicsFit_TS( object, k=k, meanThres=meanThres, s=s, d=d )            
            },
            IO = {
                # two-sample analysis (Input only)
                
                message( "Info: two-sample analysis (Input only)." )
                fit <- .mosaicsFit_IO( object, k=k, d=d, truncProb=truncProb )    
            }
        )
        
        message( "Info: done!" )
        
        return(fit)
    }
)

# MOSAiCS one-sample analysis

.mosaicsFit_OS <- function( binData, k=3, meanThres=0 )
{    
    Y <- binData@tagCount
    M <- binData@mappability
    GC <- binData@gcContent
    
    message( "Info: use adaptive griding." )
    message( "Info: fitting background model..." )    
    fitParam <- .adapGridMosaicsZ0_OS(
        Y=Y, M=M, GC=GC, min_n_MGC=50, grids_MGC=c(0.01,0.02,0.04,0.10,0.20,0.50) )
    fitZ0 <- .rlmFit_OS( parEst=fitParam, mean_thres=meanThres, Y=Y, M=M, GC=GC )
    message( "Info: done!" )
    
    message( "Info: fitting one-signal-component model..." )
    fitZ1_1S <- .mosaicsZ1_1S( fitZ0, Y=Y, k=k )
    message( "Info: fitting two-signal-component model..." )
    fitZ1_2S <- .mosaicsZ1_2S( fitZ0, Y=Y, k=k )
    
    message( "Info: calculating BIC of fitted models..." )
    fitBIC_1S <- .calcModelBIC( fitZ1=fitZ1_1S, Y=Y, k=k, model="1S", type="BIC", npar=9 )
    fitBIC_2S <- .calcModelBIC( fitZ1=fitZ1_2S, Y=Y, k=k, model="2S", type="BIC", npar=12 )
    
    mosaicsEst <- new( "MosaicsFitEst",
        pi0=fitZ0$pi0, a=fitZ0$a, betaEst=fitZ0$betaEst, muEst=fitZ0$muEst,
        b=fitZ1_1S$b, c=fitZ1_1S$c,
        p1=fitZ1_2S$p1, b1=fitZ1_2S$b1, c1=fitZ1_2S$c1, b2=fitZ1_2S$b2, c2=fitZ1_2S$c2,
        analysisType="OS" )
        
    mosaicsParam <- new( "MosaicsFitParam", k=k, meanThres=meanThres )
    
    new( "MosaicsFit",
        mosaicsEst=mosaicsEst, mosaicsParam=mosaicsParam,
        coord=binData@coord, tagCount=binData@tagCount, 
        mappability=binData@mappability, gcContent=binData@gcContent,
        bic1S=fitBIC_1S, bic2S=fitBIC_2S )
}

# MOSAiCS two-sample analysis (with M & GC)

.mosaicsFit_TS <- function( binData, k=3, meanThres=1, s=2, d=0.25 )
{
    Y <- binData@tagCount
    M <- binData@mappability
    GC <- binData@gcContent
    X <- binData@input
    
    message( "Info: use adaptive griding." )
    message( "Info: fitting background model..." ) 
    fitParam <- .adapGridMosaicsZ0_TS( Y=Y, M=M, GC=GC, X=X,
        min_n_MGC=50, grids_MGC=c(0.01,0.02,0.04,0.10,0.20,0.50), min_n_X=200 )
    fitZ0 <- .rlmFit_TS( parEst=fitParam, mean_thres=meanThres, s=s, d=d,
        Y=Y, M=M, GC=GC, X=X )
    message( "Info: done!" )
    
    message( "Info: fitting one-signal-component model..." )
    fitZ1_1S <- .mosaicsZ1_1S( fitZ0, Y=Y, k=k )
    message( "Info: fitting two-signal-component model..." )
    fitZ1_2S <- .mosaicsZ1_2S( fitZ0, Y=Y, k=k )
    
    message( "Info: calculating BIC of fitted models..." )
    fitBIC_1S <- .calcModelBIC( fitZ1=fitZ1_1S, Y=Y, k=k, model="1S", type="BIC", npar=11 )
    fitBIC_2S <- .calcModelBIC( fitZ1=fitZ1_2S, Y=Y, k=k, model="2S", type="BIC", npar=14 )
    
    mosaicsEst <- new( "MosaicsFitEst",
        pi0=fitZ0$pi0, a=fitZ0$a, betaEst=fitZ0$betaEst, muEst=fitZ0$muEst,
        b=fitZ1_1S$b, c=fitZ1_1S$c,
        p1=fitZ1_2S$p1, b1=fitZ1_2S$b1, c1=fitZ1_2S$c1, b2=fitZ1_2S$b2, c2=fitZ1_2S$c2,
        analysisType="TS" )
        
    mosaicsParam <- new( "MosaicsFitParam", k=k, meanThres=meanThres, s=s, d=d )
    
    new( "MosaicsFit",
        mosaicsEst=mosaicsEst, mosaicsParam=mosaicsParam,
        coord=binData@coord, tagCount=binData@tagCount, input=binData@input, 
        mappability=binData@mappability, gcContent=binData@gcContent,
        bic1S=fitBIC_1S, bic2S=fitBIC_2S )
}

# MOSAiCS two-sample analysis (Input only)

.mosaicsFit_IO <- function( binData, k=3, d=0.25, truncProb=0.9999 )
{
    Y <- binData@tagCount
    X <- binData@input
    
    message( "Info: fitting background model..." ) 
    fitParam <- .mosaicsZ0( Y=Y, analysisType="IO", X=X, truncProb=truncProb )
    fitZ0 <- .rlmFit_IO( parEst=fitParam, d=d, Y=Y, X=X )
    message( "Info: done!" )
    
    message( "Info: fitting one-signal-component model..." )
    fitZ1_1S <- .mosaicsZ1_1S( fitZ0, Y=Y, k=k )
    message( "Info: fitting two-signal-component model..." )
    fitZ1_2S <- .mosaicsZ1_2S( fitZ0, Y=Y, k=k )
    
    message( "Info: calculating BIC of fitted models..." )
    fitBIC_1S <- .calcModelBIC( fitZ1=fitZ1_1S, Y=Y, k=k, model="1S", type="BIC", npar=6 )
    fitBIC_2S <- .calcModelBIC( fitZ1=fitZ1_2S, Y=Y, k=k, model="2S", type="BIC", npar=9 )
    
    mosaicsEst <- new( "MosaicsFitEst",
        pi0=fitZ0$pi0, a=fitZ0$a, betaEst=fitZ0$betaEst, muEst=fitZ0$muEst,
        b=fitZ1_1S$b, c=fitZ1_1S$c,
        p1=fitZ1_2S$p1, b1=fitZ1_2S$b1, c1=fitZ1_2S$c1, b2=fitZ1_2S$b2, c2=fitZ1_2S$c2,
        analysisType="IO" )
        
    mosaicsParam <- new( "MosaicsFitParam", k=k, d=d )
    
    new( "MosaicsFit",
        mosaicsEst=mosaicsEst, mosaicsParam=mosaicsParam,
        coord=binData@coord, tagCount=binData@tagCount, input=binData@input,
        bic1S=fitBIC_1S, bic2S=fitBIC_2S )
}
