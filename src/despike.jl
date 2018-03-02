"""
    Despike{R}(reference::Function,n,k)
    (d::Despike)(q::Quantity)

Remove spikes from a Quantity

Despike applies a function (`reference`)
to a moving window of size `k` to compute
a reference time series. 

Spikes are defined as those points lying `n` times
the standard deviation of the residuals away from 
the reference time series.

What it does with the spikes depends on the 
value of the type parameter `R`

`R` can be one of 

- `:filter`: Replace the entire time series with the reference series
- `:reference`: Replace identified spikes with the reference series
- `:nan`: Replace spikes with NaN
"""
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

function threshold{T<:Quantity}(q::T,t,z)
    newQ = quantity(q)[:]
    newQ[newQ.>t] = z
    T(times(q),newQ)
end
