"""
Discharge estimates are biased, and this bias corrupts
material budgets based on them. We adjust for this 
bias by assuming that the water balance is exactly satisfied.
"""
function adjustbias(Q::Discharge,T::Quantity,M::Mask)
    F = T*Q # Observed instantaneous flux
    t,Fo = tidalbalance(F,M) # Observed material balance
    _,To = tidalbalance(T,M) # Tidal integral of T balance
    _,μ = tidalaverage(Q,M) # Observed average water balance
    Fo.-μ.*To #
end
