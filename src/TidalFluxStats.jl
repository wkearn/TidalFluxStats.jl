module TidalFluxStats

using TidalFluxQuantities

export tidalbalance,
    tidalaverage,
    segment

include("averages.jl")

end # module
