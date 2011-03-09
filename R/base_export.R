
# export all possible info

.exportTXT <- function( peakList, fileLoc, fileName, chrID="chrNA" )
{
    outFormat <- as.matrix( data.frame( chrID, peakList ) )
    
    currentDir <- getwd()
    setwd( fileLoc )
    cat( file=fileName )
    
    # variable names
    
    cat( colnames(outFormat), file=fileName, sep="\t", append=TRUE )
    cat( "\n", file=fileName, append=TRUE )
    
    # peak list
     
    for ( i in 1:nrow(outFormat) )
    {
        cat( outFormat[i,], file=fileName, sep="\t", append=TRUE )
        cat( "\n", file=fileName, append=TRUE )
    }
    setwd( currentDir )
    
    message( "Info: peak file was exported in TXT format:" )
    message( "Info: file name = ", fileName )
    message( "Info: directory = ", fileLoc )
}

# GFF format (score = peak count)

.exportGFF <- function( peakList, fileLoc, fileName, chrID="chrNA" )
{
    # GFF: seqname, source, feature, start, end, score, strand, frame, group
    
    outFormat <- as.matrix( data.frame( chrID, "MOSAiCS", "MOSAiCS_peak",
        peakList$peakStart, peakList$peakStop, peakList$aveChipCount,
        ".", ".", "." ) )
    
    currentDir <- getwd()
    setwd( fileLoc )
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\"'
    cat( line0, "\n", file=fileName )
    for ( i in 1:nrow(outFormat) )
    {
        cat( outFormat[i,], file=fileName, sep="\t", append=TRUE )
        cat( "\n", file=fileName, append=TRUE )
    }
    setwd( currentDir )
    
    message( "Info: peak file was exported in GFF format:" )
    message( "Info: file name = ", fileName )
    message( "Info: directory = ", fileLoc )
}

# BED format (score = peak count)

.exportBED <- function( peakList, fileLoc, fileName, chrID="chrNA" )
{
    # BED: (required) chrom, chromStart, chromEnd
    # BED: (optional) name, score, strand, thickStart, thinkEnd, itemRgb,
    #                   blockCount, blockSizes, blockStarts
    
    outFormat <- as.matrix( data.frame( chrID,
        peakList$peakStart, peakList$peakStop, "MOSAiCS_peak", peakList$aveChipCount ) )
    
    currentDir <- getwd()
    setwd( fileLoc )
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\" useScore=1'
    cat( line0, "\n", file=fileName )
    for ( i in 1:nrow(outFormat) )
    {
        cat( outFormat[i,], file=fileName, sep="\t", append=TRUE )
        cat( "\n", file=fileName, append=TRUE )
    }
    setwd( currentDir )
    
    message( "Info: peak file was exported in BED format:" )
    message( "Info: file name = ", fileName )
    message( "Info: directory = ", fileLoc )
}
