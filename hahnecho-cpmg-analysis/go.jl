using LinearAlgebra
using Measurements
using Plots
using Formatting
using LsqFit
using Printf


## initialisation

zscore(x::Measurement) = Measurements.value(x) / Measurements.uncertainty(x)

include("constants.jl")
include("types.jl")
include("data.jl")
include("residuals.jl")
include("plotting.jl")


## load data

# import list of methyl names - numbering here will be used to select methyls for further analysis
methyls = loadmethyls("data/methyls.txt")
nmethyls = length(methyls)
for (i, name) in zip(1:nmethyls, methyls)
    println("$i\t$name")
end

# import Hahn echo relaxation rates and specify associated field strengths
R2data_raw = loadR2("data/R2data.txt", [600, 700, 800, 950])

# import S2tc and CSA measurements
S2tc = loadS2tc("data/S2tc.txt")
ΔσH = load1Hcsa("data/etaH.txt", S2tc, 950)
ΔσC = load13Ccsa("data/13Ccsa.txt")

# subtract CSA contributions from measured R2 rates
R2data = subtractCSA(R2data_raw, S2tc, ΔσH, ΔσC)

# import CPMG data - file format is based on parsed ChemEx output
# specify experiment type, field strength, relaxation time (0.04 s) and a minimum error (0.3 s-1)
CPMGdata1 = loadCPMG("data/950-MQ-CPMG.exp", methyls, :MQ, 950, 0.04, min_err=0.3)
CPMGdata2 = loadCPMG("data/800-MQ-CPMG.exp", methyls, :MQ, 800, 0.04, min_err=0.3)
CPMGdata3 = loadCPMG("data/800-1H-CPMG.exp", methyls, :SQ1H, 800, 0.04, min_err=0.3)

## select a group of residues for analysis

# # FIRST GROUP
# active_methyls = [4, 6, 13, 34, 35, 3, 7, 33]
# spinpar0 = [
#     # ΔδH, ΔδC, R20(ZQ, DQ, QQ⁻, QQ⁺), R20(CPMG1, CPMG2, CPMG3)    
#     -0.015  0.5  5.0  5.0  50.0  50.0  10.0  10.0  10.0
#     0.015  0.4  5.0  5.0  50.0  50.0  10.0  10.0  10.0
#     0.03  0.4  5.0  5.0  50.0  50.0  10.0  10.0  10.0
#     -0.03  0.4  5.0  5.0  50.0  50.0  10.0  10.0  10.0
#     -0.02  0.3  5.0  5.0  50.0  50.0  10.0  10.0  10.0
#     .01 .1 5.0  5.0  50.0  50.0  10.0  10.0  10.0
#     .01 .1 5.0  5.0  50.0  50.0  10.0  10.0  10.0
#     .01 .1 5.0  5.0  50.0  50.0  10.0  10.0  10.0
# ]
# mp0 = modelpars_2state(900, .03)  # initial estimate for kex, pB


# SECOND GROUP
active_methyls = [8, 20, 21, 19, 22]
spinpar0 = [
    # ΔδH, ΔδC, R20(ZQ, DQ, QQ⁻, QQ⁺), R20(CPMG1, CPMG2, CPMG3)
    -0.0754704  1.15793   29.4364   12.1437   156.053   30.1853   6.66643  10.1832   12.3749
    0.0518553  0.572343   6.43592   6.26245   26.4896  52.4839   7.10874   6.46275  12.9216
   -0.1         0.247   7.6476   10.7173    40.3527  40.0117  10.4753    9.84693  16.6869
   0.0518553  0.572343   6.43592   6.26245   26.4896  52.4839   7.10874   6.46275  12.9216
   -0.0855321  0.446828   7.6476   10.7173    95.3527  73.0117  10.4753    9.84693  16.6869
    ]
mp0 = modelpars_2state(6900, .16)  # initial estimate for kex, pB


# apply selection to imported data
methyls = methyls[active_methyls]
R2data = R2data[active_methyls]
CPMGdata1 = CPMGdata1[active_methyls]
CPMGdata2 = CPMGdata2[active_methyls]
CPMGdata3 = CPMGdata3[active_methyls]
nmethyls = length(active_methyls)


## prepare functions for fitting

# get list of field strengths for HE measurements
R2_B0s = R2data[1].B0

# pack model parameters and spin parameters into a single vector
function pack(mp, spinpar)
    vec([mp.kex; mp.pB; vec(spinpar)])
end

# unpack parameter vector into model and spin parameters
function unpack(p) 
    mp = modelpars_2state(p[1:2]...)
    spinpars = reshape(p[3:end], nmethyls, :)
    return mp, spinpars
end

# simulate collection of experiments - this may need to be edited, e.g. if using more/less CPMG experiments
function simulate_spins(mp, spinpars)
    nmethyls, ~ = size(spinpars)
    simR2 = []
    simCPMG1 = []
    simCPMG2 = []
    simCPMG3 = []
    for i=1:nmethyls
        ΔδH = spinpars[i, 1]
        ΔδC = spinpars[i, 2]
        R2_ZQ = spinpars[i, 3]
        R2_DQ = spinpars[i, 4]
        R2_QQ⁻ = spinpars[i, 5]
        R2_QQ⁺ = spinpars[i, 6]
        R2cpmg1 = spinpars[i, 7]
        R2cpmg2 = spinpars[i, 8]
        R2cpmg3 = spinpars[i, 9]
        push!(simR2, sim_R2s(mp, ΔδH, ΔδC, R2_ZQ, R2_DQ, R2_QQ⁻, R2_QQ⁺, R2_B0s))
        push!(simCPMG1, sim_CPMG(CPMGdata1[i], mp, ΔδH, ΔδC, R2cpmg1))
        push!(simCPMG2, sim_CPMG(CPMGdata2[i], mp, ΔδH, ΔδC, R2cpmg2))
        push!(simCPMG3, sim_CPMG(CPMGdata3[i], mp, ΔδH, ΔδC, R2cpmg3))
    end
    return simR2, simCPMG1, simCPMG2, simCPMG3
end

# calculate residuals - this may need to be edited, e.g. if using more/less CPMG experiments
function calcresiduals(mp, spinpars)
    nmethyls, ~ = size(spinpars)

    resid = Vector{Float64}()

    for i=1:nmethyls
        ΔδH = spinpars[i, 1]
        ΔδC = spinpars[i, 2]
        R2_ZQ = spinpars[i, 3]
        R2_DQ = spinpars[i, 4]
        R2_QQ⁻ = spinpars[i, 5]
        R2_QQ⁺ = spinpars[i, 6]
        R2cpmg1 = spinpars[i, 7]
        R2cpmg2 = spinpars[i, 8]
        R2cpmg3 = spinpars[i, 9]
        append!(resid, residuals_R2(R2data[i], mp, ΔδH, ΔδC, R2_ZQ, R2_DQ, R2_QQ⁻, R2_QQ⁺, R2data[i].B0))
        append!(resid, residuals_CPMG(CPMGdata1[i], mp, ΔδH, ΔδC, R2cpmg1))
        append!(resid, residuals_CPMG(CPMGdata2[i], mp, ΔδH, ΔδC, R2cpmg2))
        append!(resid, residuals_CPMG(CPMGdata3[i], mp, ΔδH, ΔδC, R2cpmg3))
    end

    return resid
end

# calculate χ² = sum of squares of residuals, weighted by uncertainty
χ2(mp, spinpar) = sum(calcresiduals(mp, spinpar).^2)


## least-squares fitting

# run the fit
fit = LsqFit.lmfit(p -> calcresiduals(unpack(p)...), pack(mp0, spinpar0), Float64[], show_trace=true)

# extract best-fit parameters
mpfit, spinparfit = unpack(fit.param)

# calculate parameter uncertainties
mpfit_err, spinparfit_err = unpack(fit.param .± stderror(fit))

# calculate best-fit χ² value
χ2(mpfit, spinparfit)

# calculate best-fit predictions and plot
simR2, simCPMG1, simCPMG2, simCPMG3 = simulate_spins(mpfit, spinparfit)
for i=1:nmethyls
    plot(plotdata(R2data[i], simR2[i]),
        plotdata([CPMGdata1[i], CPMGdata2[i], CPMGdata3[i]],
            [simCPMG1[i], simCPMG2[i], simCPMG3[i]]),
        layout=(1,2),title=methyls[i],size=(800,300)) |> display
end

println()
println("kex = $(mpfit_err.kex) s⁻¹")
println("pB = $(mpfit_err.pB)")
spinparfit_err


## plot a single residue fit (used for main text figures)
font = Plots.font("Helvetica", 6)
plotstyle = Dict(
    :guidefont=>font,
    :xtickfont=>font,
    :ytickfont=>font,
    :legendfont=>font,
    :framestyle=>:box,
    :grid=>false,
    :lw=>0.5)

i = 1  # select residue number
plot(plotdata(R2data[i], simR2[i]),
    plotdata([CPMGdata1[i], CPMGdata2[i], CPMGdata3[i]], [simCPMG1[i], simCPMG2[i], simCPMG3[i]]),
    layout=(1,2),
    legend=:topright,
    size=(500,200);
    plotstyle...)


## combine plots (used for SI figures)
font = Plots.font("Helvetica", 6)
plotstyle = Dict(
    :guidefont=>font,
    :xtickfont=>font,
    :ytickfont=>font,
    :legendfont=>font,
    :framestyle=>:box,
    :grid=>false,
    :lw=>0.5)

plts = []
# order = [2, 7, 3, 8, 4, 6, 1, 5]  # first group
order = [1, 4, 2, 3, 5]  # second group
for j=1:nmethyls
    i = order[j]
    push!(plts, plotdata(R2data[i], simR2[i]))
    push!(plts, plotdata([CPMGdata1[i], CPMGdata2[i], CPMGdata3[i]],
            [simCPMG1[i], simCPMG2[i], simCPMG3[i]]))
end
for i=1:6 # pad second plots
    push!(plts, plot())
end
plot(plts..., layout=(4,4), legend=nothing, size=(800,600), ms=3; plotstyle...)
