
# generic methods for "BinData" class

setGeneric( "coord",
    function( object, ... )
    standardGeneric("coord")
)

setGeneric( "tagCount",
    function( object, ... )
    standardGeneric("tagCount")
)

setGeneric( "mappability",
    function( object, ... )
    standardGeneric("mappability")
)

setGeneric( "gcContent",
    function( object, ... )
    standardGeneric("gcContent")
)

setGeneric( "input",
    function( object, ... )
    standardGeneric("input")
)

# generic methods for "MosaicsFit" class

setGeneric( "mosaicsFit",
    function( object, ... )
    standardGeneric("mosaicsFit")
)

setGeneric( "estimates",
    function( object, ... )
    standardGeneric("estimates")
)

# generic methods for "MosaicsPeak" class

setGeneric( "mosaicsPeak",
    function( object, ... )
    standardGeneric("mosaicsPeak")
)

setGeneric( "peakList",
    function( object, ... )
    standardGeneric("peakList")
)

setGeneric( "export",
    function( object, ... )
    standardGeneric("export")
)

setGeneric( "bdBin",
    function( object, ... )
    standardGeneric("bdBin")
)

setGeneric( "empFDR",
    function( object, ... )
    standardGeneric("empFDR")
)
