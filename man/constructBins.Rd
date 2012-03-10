\name{constructBins}
\alias{constructBins}
\title{
Construct bin-level ChIP-sep data from an aligned read file
}
\description{
Preprocess and construct bin-level ChIP-sep data from an aligned read file.
}
\usage{
constructBins( infile=NULL, fileFormat=NULL, outfileLoc="./", 
    byChr=FALSE, excludeChr=NULL, 
    fragLen=200, binSize=200, capping=0, perl = "perl" )
}
\arguments{
  \item{infile}{
    Name of the aligned read file to be processed.
}
  \item{fileFormat}{
    Format of the aligned read file to be processed.
    Currently, \code{constructBins} permits the following aligned read file formats:
        \code{"eland_result"} (Eland result), \code{"eland_extended"} (Eland extended),
        \code{"eland_export"} (Eland export), \code{"bowtie"} (default Bowtie),
        \code{"sam"} (SAM), and \code{"bed"} (BED).
}
  \item{outfileLoc}{
    Directory of processed bin-level files. 
    By default, processed bin-level files are exported to the current directory. 
}
  \item{byChr}{
    Construct separate bin-level file for each chromosome? 
    Possible values are \code{TRUE} or \code{FALSE}.
    If \code{byChr=FALSE}, all chromosomes are exported into one file. Default is \code{FALSE}.
}
  \item{excludeChr}{
    Vector of chromosomes that are excluded from the analysis.
}  
  \item{fragLen}{
    Average fragment length. Default is 200.
}
  \item{binSize}{
    Size of bins. Default is 200.
}
  \item{capping}{
    Maximum number of reads allowed to start at each nucleotide position. 
    To avoid potential PCR amplification artifacts, the maximum number of reads
    that can start at a nucleotide position is capped at \code{capping}. 
    Capping is not applied if non-positive \code{capping} is used.
    Default is 0 (no capping).
}
  \item{perl}{
    Name of the perl executable to be called. Default is \code{"perl"}.
}
}
\details{
Bin-level files are constructred from the aligned read file and exported to \code{outfileLoc}.
If \code{byChr=FALSE}, bin-level files are named 
as \code{[infileName]_fragL[fragLen]_bin[binSize].txt},
If \code{byChr=TRUE}, bin-level files are named 
as \code{[infileName]_fragL[fragLen]_bin[binSize]_[chrID].txt},
where \code{[chrID]} is chromosome ID that reads align to.
These chromosome IDs are extracted from the aligned read file.
Constructed bin-level files can be loaded into the R environment using the method \code{readBins}.

\code{constructBins} currently supports the following aligned read file formats: 
Eland result (\code{"eland_result"}), Eland extended (\code{"eland_extended"}),
Eland export (\code{"eland_export"}), default Bowtie (\code{"bowtie"}), 
SAM (\code{"sam"}), and BED (\code{"bed"}).
This method assumes that these aligned read files are obtained from single-end tag (SET) experiments
and retains only reads mapping uniquely to the reference genome. 
}
\value{
Processed bin-level files are exported to the directory specified in \code{outfileLoc}.
}
\references{
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{readBins}}, \code{\linkS4class{BinData}}.
}
\examples{
\dontrun{
constructBins( infile="/scratch/eland/STAT1_eland_results.txt", 
    fileFormat="eland_result", outfileLoc="/scratch/eland/", 
    byChr=FALSE, excludeChr="chrM", fragLen=200, binSize=200, capping=0 )
}
}
\keyword{models}
\keyword{methods}
