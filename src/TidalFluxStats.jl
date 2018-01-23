module TidalFluxStats

using TidalFluxQuantities

export tidalbalance,
    tidalaverage,
    segment,
    superimpose,
    adjustbias

include("averages.jl")
include("superimpose.jl")
include("bias.jl")

end # module
