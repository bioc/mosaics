
.reportSummary <- function( summaryFileName, resultList,
    chipDir=NULL, chipFileName=NULL, chipFileFormat=NULL, 
    controlDir=NULL, controlFileName=NULL, controlFileFormat=NULL, 
    binfileDir=NULL, peakDir=NULL, peakFileName=NULL, peakFileFormat=NULL, 
    byChr=FALSE, FDR=0.05, fragLen=200, binSize=fragLen, capping=0, 
    analysisType="IO", d=0.25, 
    signalModel="BIC", maxgap=fragLen, minsize=50, thres=10 ) {
    
    cat( "MOSAiCS: Summary of model fitting and peak calling\n", file=summaryFileName )
    cat( "\n", file=summaryFileName, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFileName, append=TRUE )
    cat( "Input/output file settings\n", file=summaryFileName, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFileName, append=TRUE )
    cat( "\n", file=summaryFileName, append=TRUE )
    
    cat( "Directory of aligned read file (ChIP):", chipDir, "\n",
        file=summaryFileName, append=TRUE )
    cat( "Name of aligned read file (ChIP): ",chipFileName,"\n", sep="",
        file=summaryFileName, append=TRUE )
    cat( "Aligned read file format (ChIP):", chipFileFormat, "\n",
        file=summaryFileName, append=TRUE )
            
    cat( "\n", file=summaryFileName, append=TRUE )
    
    cat( "Directory of aligned read file (control):", controlDir, "\n",
        file=summaryFileName, append=TRUE )
    cat( "Name of aligned read file (control): ",controlFileName,"\n", sep="",
        file=summaryFileName, append=TRUE )
    cat( "Aligned read file format (control):", controlFileFormat, "\n",
        file=summaryFileName, append=TRUE )
                 
    cat( "\n", file=summaryFileName, append=TRUE )
    
    for ( ff in 1:length(peakFileFormat) ) {
        if ( length(peakDir) == 1 ) {
            peakDirFF <- peakDir
        } else {
            peakDirFF <- peakDir[ff]
        }
        
        cat( "Directory of peak result file:", peakDirFF, "\n",
            file=summaryFileName, append=TRUE )
        cat( "Name of peak result file: ",peakFileName[ff],"\n", sep="",
            file=summaryFileName, append=TRUE )
        cat( "Peak result file format:", peakFileFormat[ff], "\n",
            file=summaryFileName, append=TRUE )
        
        cat( "\n", file=summaryFileName, append=TRUE )
    }
               
    #cat( "\n", file=summaryFileName, append=TRUE )    
    
    cat( "------------------------------------------------------------\n", 
        file=summaryFileName, append=TRUE )
    cat( "Parameter settings\n", file=summaryFileName, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFileName, append=TRUE )
    cat( "\n", file=summaryFileName, append=TRUE )
    
    if ( byChr ) {    
        cat( "Genome-wide or chromosome-wise analysis? Chromosome-wise analysis\n",
            file=summaryFileName, append=TRUE )
    } else  {    
        cat( "Genome-wide or chromosome-wise analysis? Genome-wide analysis\n",
            file=summaryFileName, append=TRUE )
    }      
    
    if ( analysisType=="OS" ) {
        analysisTypeOut <- "One-sample analysis"
    } else if ( analysisType=="TS" ) {
        analysisTypeOut <- "Two-sample analysis (with mappability & GC content)"
    } else if ( analysisType=="IO" ) {
        analysisTypeOut <- "Two-sample analysis (Input only)"
    }
    
    if ( signalModel=="BIC" ) {
        signalModelOut <- "Automatic signal model selection based on BIC"
    } else if ( signalModel=="1S" ) {
        signalModelOut <- "One-signal-component model"
    } else if ( signalModel=="2S" ) {
        signalModelOut <- "Two-signal-component model"
    }
      
    cat( "False discovery rate (FDR):", FDR, "\n", file=summaryFileName, append=TRUE )
    cat( "Fragment length:", fragLen, "\n", file=summaryFileName, append=TRUE )
    cat( "Bin size:", binSize, "\n", file=summaryFileName, append=TRUE )
    if ( capping > 0 ) {
        cat( "Maximum number of reads allowed in each nucleotide:", capping, "\n",
            file=summaryFileName, append=TRUE )
    }
    cat( "Analysis type:", analysisTypeOut, "\n", file=summaryFileName, append=TRUE )
    cat( "d:", d, "\n", file=summaryFileName, append=TRUE )
    cat( "Signal model:", signalModelOut, "\n", file=summaryFileName, append=TRUE )
    cat( "maxgap:", maxgap, "\n", file=summaryFileName, append=TRUE )
    cat( "minsize:", minsize, "\n", file=summaryFileName, append=TRUE )
    cat( "thres:", thres, "\n", file=summaryFileName, append=TRUE )
       
    cat( "\n", file=summaryFileName, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFileName, append=TRUE )
    cat( "Peak calling summary\n", file=summaryFileName, append=TRUE )
    cat( "------------------------------------------------------------\n", 
        file=summaryFileName, append=TRUE )
    cat( "\n", file=summaryFileName, append=TRUE )
    
    outFormat <- data.frame( 
        resultList$chrID, resultList$n_peaks, 
        resultList$peak_width, resultList$opt_sig_model,
        stringsAsFactors=FALSE )
    colnames(outFormat) <- c( "chrID", "# peaks", 
        "Median peak width", "Optimal/specified signal model" )
      
    cat( as.character(colnames(outFormat)), file=summaryFileName, sep="\t", append=TRUE )
    cat( "\n", file=summaryFileName, append=TRUE )
    cat( rep("-----",3), file=summaryFileName, sep="\t", append=TRUE )
    cat( "\t\t\t-----", file=summaryFileName, sep="\t", append=TRUE )
    cat( "\n", file=summaryFileName, append=TRUE )
    
    # peak list
     
    for ( i in 1:nrow(outFormat) )
    {
        cat( as.character(outFormat[i,])[1:3], file=summaryFileName, sep="\t", append=TRUE )
        cat( "\t\t", as.character(outFormat[i,])[4], 
            file=summaryFileName, sep="\t", append=TRUE )
        cat( "\n", file=summaryFileName, append=TRUE )
    }
}
