\name{mosaicsRunAll}
\alias{mosaicsRunAll}
\title{
Analyze ChIP-seq data using the MOSAiCS framework
}
\description{
Construct bin-level ChIP-sep data from aligned read files of ChIP and control samples,
fit MOSAiCS model, call peaks, and export peak calling results and reports for diagnostics.
}
\usage{
mosaicsRunAll( chipDir=NULL, chipFileName=NULL, chipFileFormat=NULL, 
    controlDir=NULL, controlFileName=NULL, controlFileFormat=NULL, 
    binfileDir=NULL, peakDir=NULL, peakFileName=NULL, peakFileFormat=NULL,
    reportSummary=FALSE, summaryDir=NULL, summaryFileName=NULL, 
    reportExploratory=FALSE, exploratoryDir=NULL, exploratoryFileName=NULL, 
    reportGOF=FALSE, gofDir=NULL, gofFileName=NULL, byChr=FALSE,
    excludeChr=NULL, FDR=0.05, fragLen=200, binSize=fragLen, capping=0, 
    analysisType="IO", bgEst=NA, d=0.25, 
    signalModel="BIC", maxgap=fragLen, minsize=50, thres=10, parallel=FALSE, nCore=8 )
}
\arguments{
  \item{chipDir}{
    Directory of the aligned read file of ChIP sample to be processed.
}
  \item{chipFileName}{
    Name of the aligned read file of ChIP sample to be processed.
}
  \item{chipFileFormat}{
    Format of the aligned read file of ChIP sample to be processed.
    Currently, \code{mosaicsRunAll} permits the following aligned read file formats:
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        and \code{"sam"} (SAM).
}
  \item{controlDir}{
    Directory of the aligned read file of control sample to be processed.
}
  \item{controlFileName}{
    Name of the aligned read file of control sample to be processed.
}
  \item{controlFileFormat}{
    Format of the aligned read file of control sample to be processed.
    Currently, \code{mosaicsRunAll} permits the following aligned read file formats:
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        and \code{"sam"} (SAM).
}
  \item{binfileDir}{
    Directory to store processed bin-level files.
}
  \item{peakDir}{
    Directory to store the peak list generated from the analysis.
}
  \item{peakFileName}{
    Name of the peak list generated from the analysis.
}
  \item{peakFileFormat}{
    Format of the peak list generated from the analysis.
    Possible values are \code{"txt"}, \code{"bed"}, and \code{"gff"}.
}
  \item{reportSummary}{
    Report the summary of model fitting and peak calling? 
    Possible values are \code{TRUE} and \code{FALSE}. Default is \code{FALSE}. 
}
  \item{summaryDir}{
    Directory to store the summary report of model fitting and peak calling.
}
  \item{summaryFileName}{
    Name of the summary report of model fitting and peak calling.
    The summary report is a text file.
}
  \item{reportExploratory}{
    Report the exploratory analysis plots? 
    Possible values are \code{TRUE} and \code{FALSE}. Default is \code{FALSE}. 
}
  \item{exploratoryDir}{
    Directory to store the exploratory analysis plots.
}
  \item{exploratoryFileName}{
    Name of the file for exploratory analysis plots.
    The exploratory analysis results are exported as PDF.
}
  \item{reportGOF}{
    Report the goodness of fit (GOF) plots? 
    Possible values are \code{TRUE} and \code{FALSE}. Default is \code{FALSE}. 
}
  \item{gofDir}{
    Directory to store the goodness of fit (GOF) plots.
}
  \item{gofFileName}{
    Name of the file for goodness of fit (GOF) plots.
    The exploratory analysis results are exported as PDF.
}
  \item{byChr}{
    Analyze ChIP-seq data for each chromosome separately or analyze it genome-wide?
    Possible values are \code{TRUE} or \code{FALSE}.
    \code{byChr=TRUE} and \code{byChr=FALSE} mean 
    chromosome-wise and genome-wide analysis, respectively.
    Default is \code{FALSE} (genome-wide analysis).
}
  \item{excludeChr}{
    Vector of chromosomes that are excluded from the analysis.
}  
  \item{FDR}{
    False discovery rate. Default is 0.05.
}  
  \item{fragLen}{
    Average fragment length. Default is 200.
}
  \item{binSize}{
    Size of bins. By default, bin size equals to \code{fragLen} (average fragment length).
}
  \item{capping}{
    Maximum number of reads allowed to start at each nucleotide position. 
    To avoid potential PCR amplification artifacts, the maximum number of reads
    that can start at a nucleotide position is capped at \code{capping}. 
    Capping is not applied if non-positive \code{capping} is used.
    Default is 0 (no capping).
}
  \item{analysisType}{Analysis type. Currently, only "IO" is supported. }    
  \item{bgEst}{Parameter to determine background estimation approach.
    Possible values are "matchLow" (estimation using bins with low tag counts) and
    "rMOM" (estimation using robust method of moment (MOM)).
    If \code{bgEst} is not specified,
            this method tries to guess its best for \code{bgEst}, 
            based on the data provided.}    
  \item{d}{Parameter for estimating background distribution.
    Default is 0.25. } 
  \item{signalModel}{Signal model.    
    Possible values are "BIC" (automatic model selection using BIC), 
    "1S" (one-signal-component model), and
    "2S" (two-signal-component model). Default is "BIC". } 
  \item{maxgap}{Initial nearby peaks are merged if the distance (in bp) 
    between them is less than \code{maxgap}. 
    By default, \code{maxgap} equals to \code{fragLen} (average fragment length). }     
  \item{minsize}{An initial peak is removed if its width is narrower than \code{minsize}. 
    Default is 50. }
  \item{thres}{A bin within initial peak is removed if its ChIP tag counts are less than \code{thres}.
    Default is 10. }
  \item{parallel}{Utilize multiple CPUs for parallel computing 
    using \code{"multicore"} package?
    Possible values are \code{TRUE} (use \code{"multicore"}) 
    or \code{FALSE} (not use \code{"multicore"}).
    Default is \code{FALSE} (not use \code{"multicore"}).}
  \item{nCore}{Number of maximum number of CPUs used for the analysis. 
    Default is 8. }
}
\details{
This method implements the work flow to analyze ChIP-seq data using the MOSAiCS framework.
It imports aligned read files of ChIP and control samples,
process them into bin-level files, fit MOSAiCS model, call peaks, and export the peak lists. 
This method is a wrapper function of \code{constructBins}, \code{readBins},
\code{mosaicsFit}, \code{mosaicsPeak}, \code{export}, and methods of classes
\code{BinData}, \code{MosaicsFit}, and \code{MosaicsPeak}.

See the vignette of the package for the illustration of the work flow
and the description of employed methods and their options.
Exploratory analysis plots and goodness of fit (GOF) plots are generated 
using the methods \code{plot} of the classes \code{BinData} and \code{MosaicsFit}, respectively.
See the help of \code{constructBins} for details of the options \code{chipFileFormat},
\code{controlFileFormat}, \code{byChr}, \code{fragLen}, \code{binSize}, and \code{capping}.
See the help of \code{readBins} for details of the option \code{excludeChr}.
See the help of \code{mosaicsFit} for details of the options \code{analysisType},
\code{bgEst}, and \code{d}.
See the help of \code{mosaicsPeak} for details of the options \code{FDR}, 
\code{signalModel}, \code{maxgap}, \code{minsize}, and \code{thres}.
See the help of \code{export} for details of the option \code{peakFileFormat}.

When the data contains multiple chromosomes, 
parallel computing can be utilized for faster preprocessing and model fitting
if \code{parallel=TRUE} and \code{multicore} package is installed.
\code{nCore} determines number of CPUs used for parallel computing.
}
\value{
Processed bin-level files are exported to the directory specified in \code{binfileDir}.
If \code{byChr=FALSE} (genome-wide analysis), 
one bin-level file is exported for each of ChIP and control samples,
where file names are \code{[chipFileName]_fragL[fragLen]_bin[binSize].txt}
and \code{[controlFileName]_fragL[fragLen]_bin[binSize].txt}, respectively.
If \code{byChr=TRUE} (chromosome-wise analysis),
bin-level files are exported for each chromosome of each of ChIP and control samples,
where file names are \code{[chrID]_[chipFileName]_fragL[fragLen]_bin[binSize].txt}
and \code{[chrID]_[controlFileName]_fragL[fragLen]_bin[binSize].txt}
(\code{[chrID]} is chromosome ID that reads align to).
The peak list generated from the analysis are exported to the directory specified in \code{peakDir}
with the file name specified in \code{peakFileName}.
If \code{reportSummary=TRUE}, the summary of model fitting and peak calling is exported 
to the directory specified in \code{summaryDir} 
with the file name specified in \code{summaryFileName} (text file).
If \code{reportExploratory=TRUE}, the exploratory analysis plots are exported 
to the directory specified in \code{exploratoryDir}
with the file name specified in \code{exploratoryFileName} (PDF file).
If \code{reportGOF=TRUE}, the goodness of fit (GOF) plots are exported 
to the directory specified in \code{gofDir} 
with the file name specified in \code{gofFileName} (PDF file).
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{constructBins}}, \code{\link{readBins}}, 
\code{\link{mosaicsFit}}, \code{\link{mosaicsPeak}}, \code{\link{export}},
\code{\linkS4class{BinData}}, \code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsPeak}}.
}
\examples{
\dontrun{
# minimal input (without any reports for diagnostics)

mosaicsRunAll( 
    chipDir="/scratch/eland/", 
    chipFileName="STAT1_eland_results.txt", 
    chipFileFormat="eland_result", 
    controlDir="/scratch/eland/", 
    controlFileName="input_eland_results.txt", 
    controlFileFormat="eland_result", 
    binfileDir="/scratch/bin/", 
    peakDir="/scratch/peak/", 
    peakFileName="STAT1_peak_list.txt", 
    peakFileFormat="txt" )
    
# generate all reports for diagnostics  
    
mosaicsRunAll( 
    chipDir="/scratch/eland/", 
    chipFileName="STAT1_eland_results.txt", 
    chipFileFormat="eland_result", 
    controlDir="/scratch/eland/", 
    controlFileName="input_eland_results.txt", 
    controlFileFormat="eland_result", 
    binfileDir="/scratch/bin/", 
    peakDir="/scratch/peak/", 
    peakFileName="STAT1_peak_list.txt", 
    peakFileFormat="txt",
    reportSummary=TRUE, 
    summaryDir="/scratch/reports/", 
    summaryFileName="mosaics_summary.txt", 
    reportExploratory=TRUE, 
    exploratoryDir="/scratch/reports/", 
    exploratoryFileName="mosaics_exploratory.pdf", 
    reportGOF=TRUE, 
    gofDir="/scratch/reports/", 
    gofFileName="mosaics_GOF.pdf",
    byChr=FALSE,
    FDR=0.05, 
    fragLen=200, 
    capping=0, 
    parallel=FALSE, 
    nCore=8 )
}
}
\keyword{models}
\keyword{methods}
