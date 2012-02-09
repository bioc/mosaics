\name{readBins}
\alias{readBins}
\title{
Import bin-level ChIP-sep data
}
\description{
Import and preprocess 
all or subset of bin-level ChIP-sep data, including ChIP data, control data, 
mappability score, GC content score, and sequence ambiguity score.
}
\usage{
readBins( type = c("chip", "M", "GC", "N"), fileName = NULL,
    excludeChr=NULL, dataType = "unique", rounding = 100, nCore=8 )
}
\arguments{
  \item{type}{
    Character vector indicating data types to be imported.
    This vector can contain \code{"chip"} (ChIP data), \code{"input"} (input data), 
    \code{"M"} (mappability score), \code{"GC"} (GC content score), and \code{"N"} (sequence ambiguity score).
    Currently, \code{readBins} permits only the following combinations:
        \code{c("chip", "input", "M", "GC", "N")}, \code{c("chip", "M", "GC", "N")},  
        \code{c("chip", "input", "N")}, and \code{c("chip", "input")}.
    Default is \code{c("chip", "M", "GC", "N")}.
}
  \item{fileName}{
    Character vector of file names, each of which matches each element of \code{type}.
    \code{type} and \code{fileName} should have the same length and
    corresponding elements in two vectors should appear in the same order.
}
  \item{excludeChr}{
    Vector of chromosomes that are excluded from the analysis.
}  
  \item{dataType}{
    How reads were processed? Possible values are 
    either \code{"unique"} (only uniquely aligned reads were retained)
    or \code{"multi"} (reads aligned to multiple locations were also retained).
}
  \item{rounding}{
    How are mappability score and GC content score rounded?
    Default is 100 and this indicates rounding of mappability score and GC content score
    to the nearest hundredth.
}
    \item{nCore}{Number of CPUs when parallel computing is utilized. }
}
\details{
Bin-level ChIP and input data can be generated from the aligned read files for your samples (e.g., files
obtained from the ELAND aligner) using the method \code{constructBins}. 
In \code{mosaics} package companion website, \url{http://www.stat.wisc.edu/~keles/Software/mosaics/}, 
we provide preprocessed mappability score, GC content score, 
and sequence ambiguity score files for diverse reference genomes. 
Please check the website and the vignette for further details.

The imported data type constraints the analysis that can be implemented.
If \code{type=c("chip", "M", "GC", "N")}, only one-sample analysis is permitted.
If \code{type=c("chip", "input")} or \code{c("chip", "input", "N")}, 
only two-sample analysis without using mappability and GC content is possible. 
For \code{type=c("chip", "input", "M", "GC", "N")},
user can do all the one- or two-sample analysis. See also help page of \code{mosaicsFit}.

When the data contains multiple chromosomes, 
parallel computing can be utilized for faster preprocessing
if \code{multicore} package is installed.
\code{nCore} determines number of CPUs used for parallel computing.
}
\value{
Construct \code{BinData} class object.
}
\references{
Kuan, PF, D Chung, G Pan, JA Thomson, R Stewart, and S Keles (2011), 
"A Statistical Framework for the Analysis of ChIP-Seq Data", 
\emph{Journal of the American Statistical Association}, Vol. 106, pp. 891-903.
}
\author{ Dongjun Chung, Pei Fen Kuan, Sunduz Keles }
\seealso{
\code{\link{constructBins}}, \code{\link{mosaicsFit}}, \code{\linkS4class{BinData}}.
}
\examples{
\dontrun{
library(mosaicsExample)
exampleBinData <- readBins( type=c("chip","input","M","GC","N"),
    fileName=c( system.file("extdata/chip_chr21.txt", package="mosaicsExample"),
    system.file("extdata/input_chr21.txt", package="mosaicsExample"),
    system.file("extdata/M_chr21.txt", package="mosaicsExample"),
    system.file("extdata/GC_chr21.txt", package="mosaicsExample"),
    system.file("extdata/N_chr21.txt", package="mosaicsExample") ) )
}
}
\keyword{models}
\keyword{methods}
