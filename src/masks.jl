"""
Mask out data above a threshold `t`
"""
function ThresholdMask(q::Quantity,t)
    Mask(times(q),Array{Bool}(quantity(q).<t))
end

@quantity MaskInt Int

MaskInt(m::Mask) = MaskInt(times(m),Array{Int}(quantity(m)))

"""
Tidally average a mask and threshold

This provides a mask for excluding tides
from computations.

- `dm`: the mask that you want to average
- `m` : the mask that segments the tides
- `t` : the percentage of good data points required
"""
function average_mask(dm::Mask,m::Mask,t)
    tm,Tm = tidalaverage(MaskInt(dm),m)
    tm,Tm.>t
end
