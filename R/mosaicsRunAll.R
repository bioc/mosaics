mosaicsRunAll <- function( 
    chipDir=NULL, chipFileName=NULL, chipFileFormat=NULL, 
    controlDir=NULL, controlFileName=NULL, controlFileFormat=NULL, 
    binfileDir=NULL, peakDir=NULL, peakFileName=NULL, peakFileFormat=NULL,
    reportSummary=FALSE, summaryDir=NULL, summaryFileName=NULL, 
    reportExploratory=FALSE, exploratoryDir=NULL, exploratoryFileName=NULL, 
    reportGOF=FALSE, gofDir=NULL, gofFileName=NULL, byChr=FALSE,
    excludeChr=NULL, FDR=0.05, fragLen=200, binSize=fragLen, capping=0, 
    analysisType="IO", bgEst=NA, d=0.25, 
    signalModel="BIC", maxgap=fragLen, minsize=50, thres=10, nCore=8 ) {
    
    # check options: input & output (required)
    
    if ( is.null(chipDir) ) { stop( "Please specify 'chipDir'!" ) }
    if ( is.null(chipFileName) ) { stop( "Please specify 'chipFileName'!" ) }
    if ( is.null(chipFileFormat) ) { stop( "Please specify 'chipFileFormat'!" ) }
    
    if ( is.null(controlDir) ) { stop( "Please specify 'controlDir'!" ) }
    if ( is.null(controlFileName) ) { stop( "Please specify 'controlFileName'!" ) }
    if ( is.null(controlFileFormat) ) { stop( "Please specify 'controlFileFormat'!" ) }
    
    if ( is.null(peakDir) ) { stop( "Please specify 'peakDir'!" ) }
    if ( is.null(peakFileName) ) { stop( "Please specify 'peakFileName'!" ) }
    if ( is.null(peakFileFormat) ) { stop( "Please specify 'peakFileFormat'!" ) }
    
    if ( is.null(binfileDir) ) { stop( "Please specify 'binfileDir'!" ) }
  
    # check options: analysisType (fixed)
    
    if ( analysisType!="IO" ) {
        message( "Info: only input-only analysis is currently supported in 'mosaicsRunAll'." )
        message( "Info: 'analysisType' is set to 'IO'." )
    }
    
    # check options: peak list
    
    if ( length(peakDir) > 1 ) {
        if ( length(peakDir) != length(peakFileName) | length(peakDir) != length(peakFileFormat) ) {
            stop( "Length of 'peakDir' should be either one or same as those of 'peakFileName' and 'peakFileFormat'!" )
        }
    }
    if ( length(peakFileName) != length(peakFileFormat) ) {
        stop( "Lengths of 'peakFileName' and 'peakFileFormat' should be same!" )
    }
    
    # check options: reports (optional)
    
    if ( reportSummary ) {
        if ( is.null(summaryDir) ) {
            stop( "Please specify 'summaryDir'!" )
        }
        if ( is.null(summaryFileName) ) {
            stop( "Please specify 'summaryFileName'!" )
        }
    }
  
    if ( reportGOF ) {
        if ( is.null(gofDir) ) {
            stop( "Please specify 'gofDir'!" )
        }
        if ( is.null(gofFileName) ) {
            stop( "Please specify 'gofFileName'!" )
        }
    }
  
    if ( reportExploratory ) {
        if ( is.null(exploratoryDir) ) {
            stop( "Please specify 'exploratoryDir'!" )
        }
        if ( is.null(exploratoryFileName) ) {
            stop( "Please specify 'exploratoryFileName'!" )
        }
    }
  
    # construction of bin-level files

    cat( "Info: constructing bin-level files...\n" )
    
    processSet <- list()
    processSet[[1]] <- c( chipDir, chipFileName, chipFileFormat )
    processSet[[2]] <- c( controlDir, controlFileName, controlFileFormat )
    
    
    if ( is.element( "multicore", installed.packages()[,1] ) ) {
        # if "multicore" package exists, utilize parallel computing with "mclapply"
        library(multicore)
        
        mclapply( processSet, function(x) {
            constructBins( 
                infileLoc = x[1], infileName = x[2], fileFormat = x[3], 
                outfileLoc = binfileDir, byChr = byChr,
                fragLen = fragLen, binSize = binSize, capping = capping )    
        }, mc.cores=nCore )
    } else {
        # otherwise, use usual "lapply"
        
        lapply( processSet, function(x) {
            constructBins( 
                infileLoc = x[1], infileName = x[2], fileFormat = x[3], 
                outfileLoc = binfileDir, byChr = byChr,
                fragLen = fragLen, binSize = binSize, capping = capping )    
        } )
    }
    
    if ( byChr ) {
        ###############################################################
        #                                                             #
        #                chromosome-wise analysis                     #
        #                                                             #
        ###############################################################
    
        # read in bin-level files                           
        
        cat( "Info: analyzing bin-level files...\n" )
        
        currentLoc <- getwd()
        
        setwd( binfileDir )
        #list_chip <- system( 
        #    paste("ls *_",chipFileName,"_fragL",fragLen,"_bin",binSize,".txt",sep=""), 
        #    intern=TRUE ) 
        #list_control <- system( 
        #    paste("ls *_",controlFileName,"_fragL",fragLen,"_bin",binSize,".txt",sep=""), 
        #    intern=TRUE )
        list_chip <- system( 
            paste("ls * | grep '",chipFileName,"_fragL",fragLen,"_bin",binSize,".txt'",sep=""), 
            intern=TRUE ) 
        list_control <- system( 
            paste("ls * | grep '",controlFileName,"_fragL",fragLen,"_bin",binSize,".txt'",sep=""), 
            intern=TRUE )
        setwd(currentLoc)
        
        # check list of chromosomes & analyze only chromosomes 
        # that bin-level files for both chip & control exist
        
        chrID_chip <- unlist( lapply( strsplit( list_chip, paste("_",chipFileName,sep="") ), 
            function(x) x[1] ) )
        chrID_control <- unlist( lapply( strsplit( list_control, paste("_",controlFileName,sep="") ), 
            function(x) x[1] ) )        
        index_chip <- which( !is.na( match( chrID_chip, chrID_control ) ) )
        index_control <- match( chrID_chip, chrID_control )
        index_list <- list()
        for ( i in 1:length(index_chip) ) {
            index_list[[i]] <- c( index_chip[i], index_control[i] )
        }
        
        # model fitting & peak calling
        # check whether rparallel is available. if so, use it.
        
        cat( "Info: fitting MOSAiCS model & call peaks...\n" )
      
        if ( length(index_chip) < nCore ) {
            nCore <- length(index_chip)
        }
        
        if ( is.element( "multicore", installed.packages()[,1] ) ) {
            # if "multicore" package exists, utilize parallel computing with "mclapply"
            library(multicore)
            
            out <- mclapply( index_list, function(x) {    
                # read in bin-level file
                
                chip_file <- list_chip[ x[1] ]
                input_file <- list_control[ x[2] ]
                
                setwd( binfileDir )                 
                bin <- readBins( type=c("chip","input"), fileName=c(chip_file,input_file) )
                setwd( currentLoc )
          
                # fit model
                
                fit <- mosaicsFit( bin, analysisType=analysisType, bgEst=bgEst, d=d )    
                
                # call peaks
                
                if ( signalModel=="BIC" ) {
                    # if not specified, use BIC
                    
                    if ( fit@bic1S < fit@bic2S ) {
                        peak <- mosaicsPeak( fit, signalModel="1S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "One-signal-component model"
                    } else {        
                        peak <- mosaicsPeak( fit, signalModel="2S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "Two-signal-component model"
                    }
                } else {                                                  
                    peak <- mosaicsPeak( fit, signalModel=signalModel, 
                        FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )
                    if ( signalModel=="1S" ) {
                        opt_sig_model <- "One-signal-component model"
                    } else {
                        opt_sig_model <- "Two-signal-component model"
                    }
                }  
          
                # keep results
                
                peakPrint <- print(peak)
                
                return( 
                    chrID=as.character( chrID_chip[ x[1] ] ), 
                    bin=bin, fit=fit, peak=peak, peakPrint=peakPrint,
                    n_peaks=nrow(peak@peakList), peak_width=median(peak@peakList$peakSize),
                    opt_sig_model=opt_sig_model )
            }, mc.cores=nCore )
        } else {
            # otherwise, use usual "lapply"
            
            out <- lapply( index_list, function(x) {    
                # read in bin-level file
                
                chip_file <- list_chip[ x[1] ]
                input_file <- list_control[ x[2] ]
                
                setwd( binfileDir )                 
                bin <- readBins( type=c("chip","input"), fileName=c(chip_file,input_file) )
                setwd( currentLoc )
          
                # fit model
                
                fit <- mosaicsFit( bin, analysisType=analysisType, bgEst=bgEst, d=d )    
                
                # call peaks
                
                if ( signalModel=="BIC" ) {
                    # if not specified, use BIC
                    
                    if ( fit@bic1S < fit@bic2S ) {
                        peak <- mosaicsPeak( fit, signalModel="1S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "One-signal-component model"
                    } else {        
                        peak <- mosaicsPeak( fit, signalModel="2S", 
                            FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                        opt_sig_model <- "Two-signal-component model"
                    }
                } else {                                                  
                    peak <- mosaicsPeak( fit, signalModel=signalModel, 
                        FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )
                    if ( signalModel=="1S" ) {
                        opt_sig_model <- "One-signal-component model"
                    } else {
                        opt_sig_model <- "Two-signal-component model"
                    }
                }  
          
                # keep results
                
                peakPrint <- print(peak)
                
                return( 
                    chrID=as.character( chrID_chip[ x[1] ] ), 
                    bin=bin, fit=fit, peak=peak, peakPrint=peakPrint,
                    n_peaks=nrow(peak@peakList), peak_width=median(peak@peakList$peakSize),
                    opt_sig_model=opt_sig_model )
            } )
        }
        
        # summarize results
        
        peakSetFinal <- c()
        for ( i in 1:length(out) ) {
            peakSetFinal <- rbind( peakSetFinal, out[[i]]$peakPrint )
        }
        
        resultList <- list()
        resultList$chrID <- resultList$n_peaks <- 
            resultList$peak_width <- resultList$opt_sig_model <- rep( NA, length(out) )
        for ( i in 1:length(out) ) {
            resultList$chrID[i] <- out[[i]]$chrID
            resultList$n_peaks[i] <- out[[i]]$n_peaks
            resultList$peak_width[i] <- out[[i]]$peak_width
            resultList$opt_sig_model[i] <- out[[i]]$opt_sig_model
        }
    } else {
        ###############################################################
        #                                                             #
        #                   genome-wide analysis                      #
        #                                                             #
        ###############################################################
        
        # read in bin-level files                           
        
        cat( "Info: analyzing bin-level files...\n" )
        
        currentLoc <- getwd()
        
        chip_file <- paste(chipFileName,"_fragL",fragLen,"_bin",binSize,".txt",sep="")
        input_file <- paste(controlFileName,"_fragL",fragLen,"_bin",binSize,".txt",sep="")
        
        # model fitting & peak calling
        # check whether rparallel is available. if so, use it.
        
        cat( "Info: fitting MOSAiCS model & call peaks...\n" )
        
        out <- list()
        
        # read in bin-level file
        
        setwd( binfileDir )                 
        out$bin <- readBins( type=c("chip","input"), fileName=c(chip_file,input_file) )
        setwd( currentLoc )
  
        # fit model
        
        out$fit <- mosaicsFit( out$bin, analysisType=analysisType, bgEst=bgEst, d=d )    
        
        # call peaks
        
        if ( signalModel=="BIC" ) {
            # if not specified, use BIC
            
            if ( out$fit@bic1S < out$fit@bic2S ) {
                out$peak <- mosaicsPeak( out$fit, signalModel="1S", 
                    FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                out$opt_sig_model <- "One-signal-component model"
            } else {        
                out$peak <- mosaicsPeak( out$fit, signalModel="2S", 
                    FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )      
                out$opt_sig_model <- "Two-signal-component model"
            }
        } else {                                                  
            out$peak <- mosaicsPeak( out$fit, signalModel=signalModel, 
                FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres )
            if ( signalModel=="1S" ) {
                out$opt_sig_model <- "One-signal-component model"
            } else {
                out$opt_sig_model <- "Two-signal-component model"
            }
        }  
  
        # keep results
        
        peakSetFinal <- print(out$peak)
        
        peakPrint <- split( peakSetFinal, peakSetFinal$chrID )        
        out$chrID <- names(peakPrint)
        out$n_peaks <- unlist( lapply( peakPrint, nrow ) )         
        out$peak_width <- unlist( lapply( peakPrint, function(x) { median(x$peakSize) } ) )
        
        resultList <- list()
        resultList$chrID <- out$chrID
        resultList$n_peaks <- out$n_peaks
        resultList$peak_width <- out$peak_width
        resultList$opt_sig_model <- rep( out$opt_sig_model, length(resultList$chrID) )
    }
  
    # write peak calling results                      
  
    cat( "Info: writing the peak list...\n" )
  
    for ( ff in 1:length(peakFileFormat) ) {
        if ( length(peakDir) == 1 ) {
            peakDirFF <- peakDir
        } else {
            peakDirFF <- peakDir[ff]
        }
        
        if ( peakFileFormat[ff] == "txt" ) {
            .exportTXT( peakList=peakSetFinal, fileLoc=peakDirFF, fileName=peakFileName[ff] )
        } else if ( peakFileFormat[ff] == "bed" ) {
            .exportBED( peakList=peakSetFinal, fileLoc=peakDirFF, fileName=peakFileName[ff] )
        } else if ( peakFileFormat[ff] == "gff" ) {
            .exportGFF( peakList=peakSetFinal, fileLoc=peakDirFF, fileName=peakFileName[ff] )
        } else {
            stop( "Inappropriate peak file format!" )  
        }
    }
  
    # report: summary
  
    cat( "Info: generating reports...\n" )
  
    if ( reportSummary ) {    
        setwd(summaryDir)
        
        .reportSummary( summaryFileName=summaryFileName, resultList=resultList,
            chipDir=chipDir, chipFileName=chipFileName, 
            chipFileFormat=chipFileFormat, 
            controlDir=controlDir, controlFileName=controlFileName, 
            controlFileFormat=controlFileFormat, 
            binfileDir=binfileDir, 
            peakDir=peakDir, peakFileName=peakFileName, peakFileFormat=peakFileFormat, 
            byChr=byChr, FDR=FDR, fragLen=fragLen, binSize=binSize, capping=capping, 
            analysisType=analysisType, d=d, 
            signalModel=signalModel, maxgap=maxgap, minsize=minsize, thres=thres )
                
        setwd(currentLoc)
    }
        
    # GOF
  
    if ( reportGOF ) {    
        setwd(gofDir)
        
        pdf(gofFileName)
        
        if ( byChr ) {
            for ( i in 1:length(out) ) {
                chrID <- out[[i]]$chrID
                fit <- out[[i]]$fit   
                
                # chrID
                
                plot( 0, 0, type="n", axes=F, ann=F )
                text( 0, 0, chrID, cex=4 )
                
                # GOF
                
                plot(fit)
            }
        } else {
            fit <- out$fit
            plot(fit)
        }        
        
        dev.off()
        
        setwd(currentLoc)
    }
        
    # exploratory analysis
  
    if ( reportExploratory ) {
        setwd(exploratoryDir)
        
        pdf(exploratoryFileName)
        
        if ( byChr ) {
            for ( i in 1:length(out) ) {
                chrID <- out[[i]]$chrID
                bin <- out[[i]]$bin
                
                # chrID
                
                plot( 0, 0, type="n", axes=F, ann=F )
                text( 0, 0, chrID, cex=4 )
                
                # exploratory plots 
                
                plot( bin )
                if ( analysisType=="IO" ) {
                    plot( bin, plotType="input" )
                }
                if ( analysisType=="OS" ) {
                    plot( bin, plotType="M" )       
                    plot( bin, plotType="GC" )            
                }
                if ( analysisType=="TS" ) {
                    plot( bin, plotType="M" )       
                    plot( bin, plotType="GC" )          
                    plot( bin, plotType="M|input" )         
                    plot( bin, plotType="GC|input" )        
                }
            }
        } else {
            bin <- out$bin
            
            # exploratory plots 
            
            plot( bin )
            if ( analysisType=="IO" ) {
                plot( bin, plotType="input" )
            }
            if ( analysisType=="OS" ) {
                plot( bin, plotType="M" )       
                plot( bin, plotType="GC" )            
            }
            if ( analysisType=="TS" ) {
                plot( bin, plotType="M" )       
                plot( bin, plotType="GC" )          
                plot( bin, plotType="M|input" )         
                plot( bin, plotType="GC|input" )        
            }
        }
        
        dev.off()
        
        setwd(currentLoc)
    }
}
