struct Despike{R}
    reference::Function
    n
    k
end

(d::Despike{:filter}){T<:Quantity}(q::T) = TidalFluxQuantities.filter_window(q,d.reference,d.k)

function (d::Despike{:reference}){T<:Quantity}(q::T)
    refQ = TidalFluxQuantities.filter_window(q,d.reference,d.k)
    res = quantity(q).-quantity(refQ)
    σ = std(res)
    idxs = find(abs.(res).>(d.n*σ))
    newQ = quantity(q)[:]
    newQ[idxs] = quantity(refQ)[idxs]
    T(times(q),newQ)    
end

function (d::Despike{:nan}){T<:Quantity}(q::T)
    refQ = TidalFluxQuantities.filter_window(q,d.reference,d.k)
    res = quantity(q).-quantity(refQ)
    σ = std(res)
    idxs = find(abs.(res).>(d.n*σ))
    newQ = quantity(q)[:]
    newQ[idxs] = NaN
    T(times(q),newQ)
end
