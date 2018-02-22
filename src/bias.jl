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
    t,Fo.-μ.*To
end

function adjustbias_average(Q::Discharge,T::Quantity,M::Mask)
    F = T*Q
    t,Fo = tidalaverage(F,M)
    _,To = tidalaverage(T,M)
    _,μ  = tidalaverage(Q,M)
    t,Fo.-μ.*To
end

function adjustbias_discharge(Q::Discharge,M::Mask)
    t,Qs = segment(Q,M)
    _,μ = tidalaverage(Q,M)
    t,Qs.-μ
end

# We want a way to concatenate discharges that have been adjusted,
# while filling (with zeros?) in outside the mask
function adjustbias_dischargefill(Q::Discharge,M::Mask)
    irise,ifall = edgeindices(M)
    t,Qs = segment(Q,M)
    _,μ = tidalaverage(Q,M)
    Qb = Qs.-μ
    B = zeros(length(Q))
    for i in 1:length(Qb)
        B[irise[i]:ifall[i]] = Qb[i]
    end
    Discharge(times(Q),B)
end
