
# export all possible info

.exportTXT <- function( peakList, filename )
{
    outFormat <- peakList
        
    # variable names
    
    cat( file=filename )
    cat( as.character(colnames(outFormat)), file=filename, sep="\t", append=TRUE )
    cat( "\n", file=filename, append=TRUE )
    
    # peak list
     
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
        cat( "\n", file=filename, append=TRUE )
    }
    
    message( "Info: peak file was exported in TXT format:" )
    message( "Info: file name = ", filename )
}

# GFF format (score = peak count)

.exportGFF <- function( peakList, filename )
{
    # GFF: seqname, source, feature, start, end, score, strand, frame, group
    
    outFormat <- data.frame( peakList$chrID, "MOSAiCS", "MOSAiCS_peak",
        peakList$peakStart, peakList$peakStop, peakList$aveChipCount,
        ".", ".", ".", stringsAsFactors=FALSE )
    
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\"'
    cat( as.character(line0), "\n", file=filename )
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
        cat( "\n", file=filename, append=TRUE )
    }
    
    message( "Info: peak file was exported in GFF format:" )
    message( "Info: file name = ", filename )
}

# BED format (score = peak count)

.exportBED <- function( peakList, filename )
{
    # BED: (required) chrom, chromStart, chromEnd
    # BED: (optional) name, score, strand, thickStart, thinkEnd, itemRgb,
    #                   blockCount, blockSizes, blockStarts
    
    outFormat <- data.frame( peakList$chrID,
        peakList$peakStart, peakList$peakStop, "MOSAiCS_peak", 
        peakList$aveChipCount, stringsAsFactors=FALSE )
    
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\" useScore=1'
    cat( as.character(line0), "\n", file=filename )
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
        cat( "\n", file=filename, append=TRUE )
    }
    
    message( "Info: peak file was exported in BED format:" )
    message( "Info: file name = ", filename )
}
