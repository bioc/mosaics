
# export all possible info

.exportTXT <- function( peakList, fileLoc, fileName )
{
    outFormat <- peakList
    
    currentDir <- getwd()
    setwd( fileLoc )
    cat( file=fileName )
    
    # variable names
    
    cat( as.character(colnames(outFormat)), file=fileName, sep="\t", append=TRUE )
    cat( "\n", file=fileName, append=TRUE )
    
    # peak list
     
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=fileName, sep="\t", append=TRUE )
        cat( "\n", file=fileName, append=TRUE )
    }
    setwd( currentDir )
    
    message( "Info: peak file was exported in TXT format:" )
    message( "Info: file name = ", fileName )
    message( "Info: directory = ", fileLoc )
}

# GFF format (score = peak count)

.exportGFF <- function( peakList, fileLoc, fileName )
{
    # GFF: seqname, source, feature, start, end, score, strand, frame, group
    
    outFormat <- data.frame( peakList$chrID, "MOSAiCS", "MOSAiCS_peak",
        peakList$peakStart, peakList$peakStop, peakList$aveChipCount,
        ".", ".", ".", stringsAsFactors=FALSE )
    
    currentDir <- getwd()
    setwd( fileLoc )
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\"'
    cat( as.character(line0), "\n", file=fileName )
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=fileName, sep="\t", append=TRUE )
        cat( "\n", file=fileName, append=TRUE )
    }
    setwd( currentDir )
    
    message( "Info: peak file was exported in GFF format:" )
    message( "Info: file name = ", fileName )
    message( "Info: directory = ", fileLoc )
}

# BED format (score = peak count)

.exportBED <- function( peakList, fileLoc, fileName )
{
    # BED: (required) chrom, chromStart, chromEnd
    # BED: (optional) name, score, strand, thickStart, thinkEnd, itemRgb,
    #                   blockCount, blockSizes, blockStarts
    
    outFormat <- data.frame( peakList$chrID,
        peakList$peakStart, peakList$peakStop, "MOSAiCS_peak", 
        peakList$aveChipCount, stringsAsFactors=FALSE )
    
    currentDir <- getwd()
    setwd( fileLoc )
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\" useScore=1'
    cat( as.character(line0), "\n", file=fileName )
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=fileName, sep="\t", append=TRUE )
        cat( "\n", file=fileName, append=TRUE )
    }
    setwd( currentDir )
    
    message( "Info: peak file was exported in BED format:" )
    message( "Info: file name = ", fileName )
    message( "Info: directory = ", fileLoc )
}
