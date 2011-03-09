
# obtained from importing bin-level data obtained from perl codes

setClass( Class="BinData",
    representation=representation(
        coord="numeric",
        tagCount="numeric",
        mappability="numeric",
        gcContent="numeric",
        input="numeric",
        dataType="character"
    )
)

# obtained from MOSAiCS Z0 & Z1 model fit

setClass( Class="MosaicsFitEst",
    representation=representation(
        pi0="numeric",
        a="numeric",
        betaEst="numeric",
        muEst="numeric",
        b="numeric",
        c="numeric",
        p1="numeric",
        b1="numeric",
        c1="numeric",
        b2="numeric",
        c2="numeric",
        analysisType="character"
    )
)

setClass( Class="MosaicsFitParam",
    representation=representation(
        k="numeric",
        meanThres="numeric",
        s="numeric",
        d="numeric"
    )
)

setClass( Class="MosaicsFit",
    representation=representation(
        mosaicsEst="MosaicsFitEst",
        mosaicsParam="MosaicsFitParam",
        coord="numeric",
        tagCount="numeric",
        mappability="numeric",
        gcContent="numeric",
        input="numeric",
        bic1S="numeric",
        bic2S="numeric"
    )
)

# obtained from final MOSAiCS peak calling

setClass( Class="MosaicsPeakParam",
    representation=representation(
        analysisType="character",
        signalModel="character",
        FDR="numeric",
        maxgap="numeric",
        minsize="numeric",
        thres="numeric"
    )
)

setClass( Class="MosaicsPeak",
    representation=representation(
        peakList="data.frame",
        peakParam="MosaicsPeakParam",
        bdBin="numeric",
        empFDR="numeric"
    )
)
