# calculate residuals for straight line fits (of R2 vs B0²)
line(B0², β, c) = @. β*B0² + c
function resid(B0², R2, β, c)
    y = Measurements.value.(R2)
    σ = Measurements.uncertainty.(R2)
    yhat = line(B0², β, c)
    return vec((y - yhat) ./ σ)
end

# define shortcuts for passing parameter vectors
line(x, p) = line(x, p[1], p[2])
resid(B0², R, p) = resid(B0², R, p[1], p[2])


# calculate χ² score using all ZQ, DQ, DQ' and QQ measurements
# performs least squares fitting to determine intercepts for each line
function calcχ2(ξH, ξC, resi)
    βZQ = Measurements.value((ξCcsa[resi]-ξHcsa[resi])^2 + (ξC-ξH)^2)
    βDQ = Measurements.value((ξCcsa[resi]+ξHcsa[resi])^2 + (ξC+ξH)^2)
    βDQQ = Measurements.value((ξCcsa[resi]-3ξHcsa[resi])^2 + (ξC-3ξH)^2)
    βQQ = Measurements.value((ξCcsa[resi]+3ξHcsa[resi])^2 + (ξC+3ξH)^2)
    
    c0 = [10.0]
    χ2 = 0.0
    fit = LsqFit.lmfit(c -> resid(B0², zq[resi,:], βZQ, c), c0, Float64[])
    χ2 += sum(fit.resid.^2)
    fit = LsqFit.lmfit(c -> resid(B0², dq[resi,:], βDQ, c), c0, Float64[])
    χ2 += sum(fit.resid.^2)
    fit = LsqFit.lmfit(c -> resid(B0², dqq[resi,:], βDQQ, c), c0, Float64[])
    χ2 += sum(fit.resid.^2)
    fit = LsqFit.lmfit(c -> resid(B0², qq[resi,:], βQQ, c), c0, Float64[])
    χ2 += sum(fit.resid.^2)
end


function fitβ(B0², R2)
    p0 = [.1, 10] # slope, intercept
    fit = LsqFit.lmfit(p -> resid(B0², R2, p), p0, Float64[])
    β = fit.param[1] ± stderror(fit)[1]
end
