module TidalFluxStats

using TidalFluxQuantities

export tidalbalance,
    tidalaverage,
    segment,
    superimpose

include("averages.jl")
include("superimpose.jl")

end # module
