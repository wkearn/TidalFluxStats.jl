# Tidal averages

"""
Take a mask and return the indices of the rising and falling edges
"""
function edgeindices(mask)
    dm = [0;diff(mask)]
    irise = find(dm.==1)
    ifall = find(dm.==-1)
    if length(irise) .!= length(ifall)
        warn("The number of upcrossings does not match the number of downcrossings")
    end
    irise,ifall
end

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

"""
Average f sampled at times ts with a trapezoid rule
"""
function average(f,ts)
    fi = integrate(f,ts)
    fi./Dates.value.(Dates.Second.(ts[end]-ts[1]))
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
    tmid,avg
end

function tidalaverage(f,mask)
    ts,fs = unzip(f)
    irise,ifall = edgeindices(mask)
    trise,tfall = edgetimes(mask,ts)
    tmid = trise + (tfall-trise)/2
    avg = [average(fs[irise[i]:ifall[i]],
                     ts[irise[i]:ifall[i]]) for i in eachindex(irise)]
    tmid,avg
end
