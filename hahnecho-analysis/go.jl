using Measurements
using Plots
using LsqFit
using ColorSchemes

val(x) = Measurements.value(x)

include("constants.jl")
include("data.jl")

## initialise data

methyls = loadmethyls("data/methyls.txt")

zq, dq, dqq, qq = loadR2("data/R2data.txt")
B0 = [600, 700, 800, 950] * 1e6 * 2π / 2.675e8;
B0² = B0.^2;

S2tc = loadS2tc("data/S2tc.txt")
ΔσH = load1Hcsa("data/etaH.txt", S2tc, 950)
ΔσC = load13Ccsa("data/13Ccsa.txt")


## calculate ξH and ξC contributions from CSA

ξHcsa = @. 2/(3*sqrt(5)) * γH * ΔσH * 1e-6 * sqrt(S2tc * 1e-9)
ξCcsa = @. 2/(3*sqrt(5)) * γC * ΔσC * 1e-6 * sqrt(S2tc * 1e-9)


## define fitting functions
# this needs to be called later as the functions rely on definitions of CSA
include("fitting.jl")
include("plotting.jl")


## select residue and run analysis
resi = 21


# get best-fit slopes (without correcting for CSA)
βZQ = fitβ(B0², zq[resi, :])
βDQ = fitβ(B0², dq[resi, :])
βDQQ = fitβ(B0², dqq[resi, :])
βQQ = fitβ(B0², qq[resi, :])

# compute χ2 surface
hval = LinRange(-.3,.3,100)
cval = LinRange(-.3,.3,120)
χ2 = zeros(length(hval), length(cval))
for j=1:length(hval)
    for k=1:length(cval)
        χ2[j,k] = calcχ2(hval[j], cval[k], resi)
    end
end

# get minimum of χ2 surface (gridsearch)
@show χ2min = minimum(vec(χ2))
@show ξHmin = hval[argmin(χ2)[1]]
@show ξCmin = cval[argmin(χ2)[2]]


# plot results
plt1 = plotlinefits(ξHmin, ξCmin, resi)
plt2 = plotχ2(hval, cval, χ2, βZQ, βDQ, βDQQ, βQQ)

plot(plt1, plt2, layout=(1,2), size=(800,400))
