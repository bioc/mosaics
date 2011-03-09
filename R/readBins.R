
# read bin-level data & process it

readBins <- function( type=c("chip","M","GC","N"), fileName=NULL, dataType='unique', rounding=100 )
{
    # [Note] Assumption: chip, input, M, GC, & N have corresponding coordinates 
    
    # error treatment
    
    .validType( type )
    .validLocation( type=type, fileName=fileName )
    .validDataType( dataType )
    
    # existence
    
    existInput <- ( length(which(type=="input")) > 0 )
    existM <- ( length(which(type=="M")) > 0 )
    existGC <- ( length(which(type=="GC")) > 0 )
    existN <- ( length(which(type=="N")) > 0 )
    
    # read bin-level data
    
    message( "Info: reading and preprocessing bin-level data..." )
    
    chipFileName <- fileName[ type=="chip" ]
    chip <- read.table( chipFileName, header=FALSE, sep='\t' )[,-1]
    
    if ( existM )
    {
        mapScoreFileName <- fileName[ type=="M" ]
        mapScore <- read.table( mapScoreFileName, header=FALSE )
    }
    if ( existGC ) {    
        gcScoreFileName <- fileName[ type=="GC" ]
        gcScore <- read.table( gcScoreFileName, header=FALSE )
    }   
    if ( existN ) {
        nNucFileName <- fileName[ type=="N" ]
        nNuc <- read.table( nNucFileName, header=FALSE )
    }    
    if ( existInput )
    {
        inputFileName <- fileName[ type=="input" ]
        input <- read.table( inputFileName, header=FALSE, sep='\t' )[,-1]
    }
    
    # error treatment
    
    if ( existInput )
    {
        if ( chipFileName == inputFileName ) {
            warning( "the same data was used for both ChIP & Input." )
            warning( "parameters cannot be properly estimated in this case." )
            warning( "remember that you will get errors in 'mosaicsFit' method!" )
        }
    } 
    
    # process bin-level data
    
    if ( existInput )
    {        
        if ( existM & existGC & existN )
        {            
            # two-sample analysis (if M, GC, N are available)
            
            outBin <- .processBin_MGCX( chip=chip, input=input,
                mapScore=mapScore, gcScore=gcScore, nNuc=nNuc, dataType=dataType, rounding=rounding )
        } else if ( existN ) {
    			  # two-sample analysis (input only, with N)
            
            outBin <- .processBin_XN( chip=chip, input=input, nNuc=nNuc, dataType=dataType )
        } else
        {
            # two-sample analysis (input only, without N)
            
            #cat( "warning: no mappability & sequence ambiguity are provided.\n" )
            #cat( "warning: processing using mappability & sequence ambiguity is skipped.\n" )
            outBin <- .processBin_X( chip=chip, input=input, dataType=dataType ) 
        }
    } else
    {
        # one-sample analysis
        
        if ( existM & existGC & existN ) {        
		        outBin <- .processBin_MGC( chip=chip,
		            mapScore=mapScore, gcScore=gcScore, nNuc=nNuc, dataType=dataType, rounding=rounding )
        } else {
        		stop( "All of mappability, GC content, and sequence ambiguity should be provided for one-sample analysis!" )	 
        }
    }
    
    message( "Info: done!\n" )
    
    # info about preprocessing
    
    if( existN ) {
		    nRatio <- length(which(nNuc[,2]==1)) / length(nNuc[,2])
		    cat( "------------------------------------------------------------\n" )
		    cat( "Info: preprocessing summary\n" )
		    cat( "------------------------------------------------------------\n" )
		    cat( "- percentage of bins with ambiguous sequences: ",round(nRatio*100),"%\n", sep="" )
		    cat( "  (these bins will be excluded from the analysis)\n" )
		    cat( "- before preprocessing:\n" )
		    cat( "\tfirst coordinates = ",min(chip[,1]),
		        ", last coordinates = ",max(chip[,1]), "\n", sep="" )
		    cat( "- after preprocessing:\n" )
		    cat( "\tfirst coordinates = ",min(outBin@coord),
		        ", last coordinates = ",max(outBin@coord), "\n", sep="" )
		    cat( "------------------------------------------------------------\n" )
    }
    
    return(outBin)
}

.validType <- function( type )
{
    # error treatment: check invalid type
    
    allType <- c("chip","M","GC","N","input")
    invalidType <- TRUE
    for ( i in 1:length(type) )
    {
        if ( length(which(!is.na(match(type[i],allType))))==0 )
        {
            invalidType <- FALSE
        }
    }
    if ( !invalidType )
    {
        stop( "Invalid 'type'! Choose among 'chip','M','GC','N','input'!" )
    }
    
    # error treatment: check whether type contains at least chip+M+GC+N or chip+input
    
    minType1 <- c("chip","M","GC","N")
    minType2 <- c("chip","input")
    
    notMin1 <- ( length(which(!is.na(match(minType1,type)))) < length(minType1) )
        # minimum requirement for one-sample analysis is not satisfied
    notMin2 <- ( length(which(!is.na(match(minType2,type)))) < length(minType2) )
        # minimum requirement for two-sample analysis (Input only) is not satisfied
    
    if ( notMin1 && notMin2 )
    {
        stop( "Minimum requirements of 'type':\nplease provide at least either c('chip','M','GC','N') or c('chip','input')!" )
    }
}

.validLocation <- function( type, fileName )
{    
    # error treatment: check whether 'fileName' info is provided
    
    if ( is.null(fileName) )
    {
        stop( "Please provide 'fileName' information!" )
    }
    
    # error treatment: check whether 'fileName' info matches 'type'    
    
    if ( length(type)!=length(fileName) )
    {
        #print(type)
        #print(fileName)
        stop( "Length of 'type' & length of 'fileName' do not match!" )
    }    
}
 
.validDataType <- function( dataType )
{   
    # error treatment: check invalid dataType
    
    if ( dataType!='unique' & dataType!='multi' )
    {
        stop( "Invalid 'dataType'! Choose either 'unique' or 'multi'!" )
    }
}

.processBin_MGC <- function( chip, mapScore, gcScore, nNuc, dataType, rounding )
{        
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(mapScore), nrow(gcScore) ) )
    dataList <- list( chip, mapScore, gcScore )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    mapScore <- mapScore[ !is.na(match(mapScore[,1],dataMinID[,1])), ]
    gcScore <- gcScore[ !is.na(match(gcScore[,1],dataMinID[,1])), ]
    nNuc <- nNuc[ !is.na(match(nNuc[,1],dataMinID[,1])), ]
    
    # exclude bins with ambiguous sequences
    
    nRegID <- which(nNuc[,2]==1)
    if ( length(nRegID)>0 ) {
	    chip <- chip[-nRegID,]
	    mapScore <- mapScore[-nRegID,]
	    gcScore <- gcScore[-nRegID,]
	    nNuc <- nNuc[-nRegID,]
    }
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    M <- round(mapScore[,2]*rounding)/rounding
    denom <- 1 - nNuc[,2]
    GC <- round(gcScore[,2]*rounding/denom)/rounding    
    
    if( dataType == 'unique' )
    {
        # in case of unique match
        if( length(which(Y>0&M==0)) > 0 )
        {
            Y[ Y>0 & M==0 ] <- 0    
        }
    }
    
    new( "BinData", coord=coord, tagCount=Y, mappability=M, gcContent=GC, dataType=dataType )
}

.processBin_X <- function( chip, input, dataType )
{            
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(input) ) )
    dataList <- list( chip, input )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    input <- input[ !is.na(match(input[,1],dataMinID[,1])), ]    
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    X <- round(input[,2])
    
    new( "BinData", coord=coord, tagCount=Y, input=X, dataType=dataType )
}

.processBin_XN <- function( chip, input, nNuc, dataType )
{            
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(input) ) )
    dataList <- list( chip, input )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    input <- input[ !is.na(match(input[,1],dataMinID[,1])), ]
    nNuc <- nNuc[ !is.na(match(nNuc[,1],dataMinID[,1])), ]
    
    # exclude bins with ambiguous sequences
    
    nRegID <- which(nNuc[,2]==1)
    if ( length(nRegID)>0 ) {
	    chip <- chip[-nRegID,]
	    input <- input[-nRegID,]
	    nNuc <- nNuc[-nRegID,]
    }
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    X <- round(input[,2]) 
    
    new( "BinData", coord=coord, tagCount=Y, input=X, dataType=dataType )
}

.processBin_MGCX <- function( chip, input, mapScore, gcScore, nNuc, dataType, rounding )
{
    # choose data with min # row
    
    minID <- which.min( c( nrow(chip), nrow(input), nrow(mapScore), nrow(gcScore) ) )
    dataList <- list( chip, input, mapScore, gcScore )
    dataMinID <- dataList[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    
    chip <- chip[ !is.na(match(chip[,1],dataMinID[,1])), ]
    input <- input[ !is.na(match(input[,1],dataMinID[,1])), ]
    mapScore <- mapScore[ !is.na(match(mapScore[,1],dataMinID[,1])), ]
    gcScore <- gcScore[ !is.na(match(gcScore[,1],dataMinID[,1])), ]
    nNuc <- nNuc[ !is.na(match(nNuc[,1],dataMinID[,1])), ]
    
    # exclude bins with ambiguous sequences
    
    nRegID <- which(nNuc[,2]==1)
    if ( length(nRegID)>0 ) {
    	chip <- chip[-nRegID,]
	    mapScore <- mapScore[-nRegID,]
	    gcScore <- gcScore[-nRegID,]
	    nNuc <- nNuc[-nRegID,]
	    input <- input[-nRegID,]
    }
    coord <- chip[,1]
    
    # round values
    
    Y <- round(chip[,2])    # round tag count in the case of multi match
    M <- round(mapScore[,2]*rounding)/rounding
    denom <- 1 - nNuc[,2]
    GC <- round(gcScore[,2]*rounding/denom)/rounding
    X <- round(input[,2])    
    
    if( dataType == 'unique' )
    {
        # in case of unique match
        if( length(which(Y>0&M==0)) > 0 )
        {
            Y[ Y>0 & M==0 ] <- 0    
        }
    }
    
    new( "BinData", coord=coord, tagCount=Y, mappability=M, gcContent=GC, input=X, dataType=dataType )
}
