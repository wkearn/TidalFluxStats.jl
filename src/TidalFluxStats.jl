module TidalFluxStats

using TidalFluxQuantities

export tidalbalance,
    tidalaverage,
    segment,
    superimpose,
    adjustbias,
    adjustbias_average,
    adjustbias_discharge,
    prism,
    adjustbias_prism,
    Despike

include("averages.jl")
include("superimpose.jl")
include("bias.jl")
include("prism.jl")
include("despike.jl")
include("masks.jl")

end # module
