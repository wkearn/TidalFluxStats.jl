"""
Scale a segmented quantity by its mean and sd
"""
function superimpose(f::Function,q::Quantity,mask::Mask)
    _,qs = segment(q,mask)
    map(f,qs)   
end

zscore(x) = (x-mean(x))/std(x)

superimpose(q::Quantity,mask::Mask) = superimpose(zscore,q,mask)

function superimpose(f::Function,q::Quantity,h::Stage,mask::Mask)
    qs = segment(q,mask)
    ts,hs = segment(h,mask)
    ms = map(x->findmax(x)[2],hs)
    maxtimes = [ts[i][ms[i]] for i in eachindex(ms)]
    ts.-maxtimes,superimpose(f,q,mask)
end
