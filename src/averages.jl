# Tidal averages

struct TidalAverage{T<:Quantity} <: Quantity{Float64}
    ts
    q
end

"""
Take a mask and return the indices of the rising and falling edges
"""
function edgeindices(mask)
    dm = [0;diff(mask)]
    irise = find(dm.==1)
    ifall = find(dm.==-1)
    # Drop the first ifall if it is before the first irise
    if ifall[1]<irise[1]
        ifall = ifall[2:end]
    end
    # Drop the last irise if it is after the last ifall
    if irise[end]>ifall[end]
        irise = irise[1:end-1]
    end
    irise,ifall
end

edgeindices(mask::Mask) = edgeindices(quantity(mask))

"""
Take a mask and return the times of the rising and falling edges
"""
function edgetimes(mask,ts)
    irise,ifall = edgeindices(mask)
    trise = ts[irise]
    tfall = ts[ifall]
    trise,tfall
end

"""
Integrate f sampled at times ts with a trapezoid rule
"""
function integrate(f,ts)
    dt = Dates.value.(Dates.Second.(diff(ts)))
    f1 = f[1:end-1]
    f2 = f[2:end]
    sum(dt.*(f1+f2)/2)
end

integrate(q::Quantity) = integrate(quantity(q),times(q))

"""
Average f sampled at times ts with a trapezoid rule
"""
function average(f,ts)
    fi = integrate(f,ts)
    fi./Dates.value.(Dates.Second.(ts[end]-ts[1]))
end

average(q::Quantity) = average(quantity(q),times(q))

function segment(q,mask)
    ts,qs = unzip(q)
    irise,ifall = edgeindices(mask)
    st = [ts[irise[i]:ifall[i]] for i in eachindex(irise)]
    sq = [qs[irise[i]:ifall[i]] for i in eachindex(irise)]
    st,sq
end

"""
Integrate a flux over tides given in the mask.

Returns midpoint times for the each tide and the accumulated flux over that tide
"""
function tidalbalance(f,mask)
    ts,fs = unzip(f)
    irise,ifall = edgeindices(mask)
    trise,tfall = edgetimes(mask,ts)
    tmid = trise + (tfall-trise)/2
    avg = [integrate(fs[irise[i]:ifall[i]],
                     ts[irise[i]:ifall[i]]) for i in eachindex(irise)]
    TidalAverage(tmid,avg)
end

function tidalaverage{T<:Quantity}(f::T,mask)
    ts,fs = unzip(f)
    irise,ifall = edgeindices(mask)
    trise,tfall = edgetimes(mask,ts)
    tmid = trise + (tfall-trise)/2
    avg = [average(fs[irise[i]:ifall[i]],
                   ts[irise[i]:ifall[i]]) for i in eachindex(irise)]
    TidalAverage{T}(tmid,avg)
end


"""
Compute the flood-ebb differential of f

q determines flood vs. ebb conditions
"""
function flood_ebb{T<:Quantity}(f::T,q,mask)
    ts,fs = unzip(f)
    f2 = fs[:]
    f2[quantity(q).<0] .*= -1
    tidalaverage(T(ts,f2),mask)
end
