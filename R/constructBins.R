
# read alignment files and construct bin-level files

constructBins <- function( infileLoc=NULL, infileName=NULL, fileFormat=NULL, 
    outfileLoc=infileLoc, byChr=FALSE, fragLen=200, binSize=fragLen, capping=0, perl = "perl" )
{   
    # preprocessing perl script embedded in "mosaics/inst/Perl/" (for SET)
    
    script <- "process_readfiles.pl"
    allFormat <- c( "eland_result", "eland_extended", "eland_export", "bowtie", "sam", "bed", "csem" )
    allFormatName <- c( "Eland result", "Eland extended", "Eland export", "Bowtie default", "SAM", "BED", "CSEM" )
    
    # Check whether perl exists
    
    CMD <- paste( perl, "-v" )
    res <- system( CMD, intern = TRUE, ignore.stderr = TRUE )
  
    if ( length(res) == 0 ) {
        # cannot proceed if perl does not exist
        
        stop( "Perl is not found on your system! Either check $PATH if installed or please install Perl." )
    } else {
    # process read files into bin-level files if perl exists
    
    # check whether minimal options are missing
    
    if ( length(infileLoc) != 1 || is.null(infileLoc) )
    {
        stop( "Please specify the location of the aligned read file!" )
    }   
    
    if ( length(infileName) != 1 || is.null(infileName) )
    {
        stop( "Please specify the name of the aligned read file!" )
    }       
    
    if ( length(fileFormat) != 1 || is.null(fileFormat) )
    {
        stop( "Please specify aligned read file format! Read '?constructBins' for supported file formats" )
    }   
    
    # check file format specification
    
    if ( length(which(!is.na(match( fileFormat, allFormat )))) == 0 )
    {
        stop( "Unsupported aligned read file format! Read '?constructBins' for supported file formats" )
    }
        
    
    # print out processing settings:
    # by default, set fragment length = 200, bin size = fragment length, capping = 3.
    
    fileFormatName <- allFormatName[ match( fileFormat, allFormat ) ]
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: setting summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "Directory of aligned read file:", infileLoc, "\n" )
    cat( "Name of aligned read file: ",infileName,"\n", sep="" )
    cat( "Aligned read file format:", fileFormatName, "\n" )
    cat( "Directory of processed bin-level files:", outfileLoc, "\n" )
    cat( "Construct bin-level files by chromosome?", ifelse(byChr,"Y","N"), "\n" )
    cat( "Fragment length:", fragLen, "\n" )
    cat( "Bin size:", binSize, "\n" )
    if ( capping > 0 ) {
        cat( "Maximum number of reads allowed in each nucleotide:", capping, "\n" )
    }
    cat( "------------------------------------------------------------\n" )
    
    
    # get path to the perl code (unified script for all file formats)
    
    Fn.Path <- system.file( file.path("Perl",script), package="mosaics")
    
    
    # process read file to bin-level files using perl codes
    
    message( "Info: reading the aligned read file and processing it into bin-level files..." )
    
    if ( capping <= 0 ) {
        capping <- 0
    }
    
    CMD <- paste( perl, 
        " ", Fn.Path,
        " ", fileFormat, 
        " ", fragLen, 
        " ", binSize, 
        " ", capping, 
        " ", infileLoc, 
        " ", outfileLoc, 
        " ", infileName, 
        " ", ifelse(byChr,"Y","N"), sep="" )
    
    res <- system( CMD, intern = TRUE )
    
    message( "Info: done!" )
    
    
    # print out processing results
    
    currentLoc <- getwd()
    setwd(outfileLoc)
    if ( byChr ) {
        #outfileName <- system( paste("ls *_",infileName,"_fragL",fragLen,"_bin",binSize,".txt",sep=""),
        #    intern=TRUE )
        outfileName <- system( paste("ls * | grep '",infileName,"_fragL",fragLen,"_bin",binSize,".txt'",sep=""),
            intern=TRUE )
    } else {
        outfileName <- paste(infileName,"_fragL",fragLen,"_bin",binSize,".txt",sep="")
    }
    setwd(currentLoc)
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: processing summary\n" )
    cat( "------------------------------------------------------------\n" )    
    cat( "Directory of processed bin-level files:", outfileLoc, "\n" )
    if ( byChr ) {
        cat( "List of processed bin-level files:\n" )
        for ( i in 1:length(outfileName) ) {                
            cat( "- ",outfileName[i],"\n", sep="" )
        }
    } else {
        cat( "Processed bin-level file: ",outfileName,"\n", sep="" )   
    }
    cat( "------------------------------------------------------------\n" )
    }
}
