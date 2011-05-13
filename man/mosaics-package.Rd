\name{mosaics-package}
\alias{mosaics-package}
\alias{mosaics}
\docType{package}
\title{
MOSAiCS (MOdel-based one and two Sample Analysis and Inference for ChIP-Seq)
}
\description{
This package provides functions for fitting MOSAiCS,
a statistical framework to analyze one-sample or two-sample ChIP-seq data.
}
\details{
\tabular{ll}{
Package: \tab mosaics\cr
Type: \tab Package\cr
Version: \tab 1.0.1\cr
Date: \tab 2011-05-05\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains three main classes, \code{BinData}, \code{MosaicsFit}, and \code{MosaicsPeak},
which represent bin-level ChIP-seq data, MOSAiCS model fit, and MOSAiCS peak calling results, respectively.
This package contains three main methods,\code{readBins}, \code{mosaicsFit}, and \code{mosaicsPeak}.
\code{readBins} method imports bin-level data and construct \code{BinData} class object.
\code{mosaicsFit} method fits MOSAiCS model using \code{BinData} class object and
 constructs \code{MosaicsFit} class object.
\code{mosaicsPeak} method calls peaks using \code{MosaicsFit} class object and
 construct \code{MosaicsPeak} class object.
\code{MosaicsPeak} class object can be exported as text files or transformed into data frame
and can be used for the downstream analysis.
This package also provides methods for simple exploratory analysis.

The \code{mosaics} package companion website, \url{http://www.stat.wisc.edu/~keles/Software/mosaics/},
provides preprocessing scripts, preprocessed files for diverse reference genomes,
and easy-to-follow instructions.
We encourage questions or requests regarding \code{mosaics} package to be posted on
our Google group, \url{http://groups.google.com/group/mosaics_user_group}.
Please check the vignette for further details on the \code{mosaics} package and these websites.
}
\author{
Dongjun Chung, Pei Fen Kuan, Sunduz Keles

Maintainer: Dongjun Chung <chungdon@stat.wisc.edu>
}
\references{
Kuan, PF, D Chung, JA Thomson, R Stewart, and S Keles (2010), "A Statistical Framework
for the Analysis of ChIP-Seq Data", submitted (\url{http://works.bepress.com/sunduz_keles/19/}).
}
\keyword{ package }
\seealso{
\code{\link{readBins}}, \code{\link{mosaicsFit}}, \code{\link{mosaicsPeak}} ,
\code{\linkS4class{BinData}}, \code{\linkS4class{MosaicsFit}}, \code{\linkS4class{MosaicsPeak}}.
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

exampleBinData
print(exampleBinData)[1:10, ]
plot(exampleBinData)
plot( exampleBinData, plotType="M" )
plot( exampleBinData, plotType="GC" )
plot( exampleBinData, plotType="input" )
plot( exampleBinData, plotType="M|input" )
plot( exampleBinData, plotType="GC|input" )

exampleFit <- mosaicsFit( exampleBinData, analysisType="TS" )

exampleFit
plot(exampleFit)
estimates(exampleFit)

examplePeak <- mosaicsPeak( exampleFit, signalModel = "2S", FDR = 0.05 )

examplePeak
print(examplePeak)[1:10, ]
export( examplePeak, type = "txt", fileLoc = "./", fileName = "TSpeakList.txt", chrID = "chr21" )
export( examplePeak, type = "bed", fileLoc = "./", fileName = "TSpeakList.bed", chrID = "chr21" )
export( examplePeak, type = "gff", fileLoc = "./", fileName = "TSpeakList.gff", chrID = "chr21" )
}
}
