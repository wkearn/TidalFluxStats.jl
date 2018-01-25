# Calculate the tidal prism
function prism(Q,ts)
    Qc = Q[:]
    Qc[Q.<0] = 0.0
    TidalFluxStats.integrate(Qc,ts)
end

function prism(Q::Discharge,M::Mask)
    st,sq = segment(Q,M)
    [prism(sq[i],st[i]) for i in eachindex(st)]
end

function adjustbias_prism(Q::Discharge,M::Mask)
    t,q = adjustbias_discharge(Q,M)
    [prism(q[i],t[i]) for i in eachindex(q)]
end
