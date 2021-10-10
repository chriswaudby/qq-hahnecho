# residual calculations


function Rex(mp::modelpars_2state, Δω)
    kex = mp.kex
    pB = mp.pB
    @. real(0.5*(kex+1im*Δω-sqrt((kex+1im*Δω)^2-4im*kex*pB*Δω)))
end

"""
simulate R2 hahn-echo relaxation rates, for 1H coherence order pH, and list of field strengths B0
"""
function sim_R2(mp::modelpars_2state, ΔδH, ΔδC, R20, pH, B0)
    Δω = 1e-6 * (ΔδC*γC + pH*ΔδH*γH) .* B0
    return R20 .+ Rex(mp, Δω)
end

function sim_R2s(mp, ΔδH, ΔδC, R2_ZQ, R2_DQ, R2_QQ⁻, R2_QQ⁺, B0s)
    dq = sim_R2(mp, ΔδH, ΔδC, R2_DQ, +1, B0s)
    zq = sim_R2(mp, ΔδH, ΔδC, R2_ZQ, -1, B0s)
    dqq = sim_R2(mp, ΔδH, ΔδC, R2_QQ⁻, -3, B0s)
    qq = sim_R2(mp, ΔδH, ΔδC, R2_QQ⁺, +3, B0s)
    return R2Data(zq, dq, dqq, qq, B0s)
end


"""
simulate CPMG relaxation rates
"""
function sim_CPMG(expt::CPMGData_MQ, mp::modelpars_2state, ΔδH, ΔδC, R20)
    # basis for MQ CPMG experiment is (ZQ⁻, ZQ⁺, DQ⁻, DQ⁺) where +/- indicates 1H coherence order
    # initial magnetisation is (1/2) (ZQ⁺ + DQ⁺)
    # final magnetisation read-out is (1/2) (ZQ⁻ + DQ⁻)

    # pulse program parameters
    T = expt.Trelax
    B0 = expt.B0
    νCPMG = expt.νCPMG
    kex = mp.kex
    pB = mp.pB

    ΔωC = γC * B0 * ΔδC * 1e-6
    ΔωH = γH * B0 * ΔδH * 1e-6
    ΔωZQ = ΔωH - ΔωC
    ΔωDQ = ΔωH + ΔωC

    # initial magnetisation
    M0 = 0.5 * [0, (1-pB), 0, (1-pB), 0, pB, 0, pB]
    Mobs = [1 0 1 0 1 0 1 0]

    # set up propagators
    UπH = [
        0 0 0 1
        0 0 1 0
        0 1 0 0
        1 0 0 0
        ] |> float
    UπC = [
        0 0 1 0
        0 0 0 1
        1 0 0 0
        0 1 0 0
        ] |> float
    UπH = kron(I(2), UπH)
    UπC = kron(I(2), UπC)
    kAB = pB*kex
    kBA = (1-pB)*kex
    K = kron([-kAB kBA
               kAB -kBA
               ], I(4))
    Ω = Diagonal([0, 0, 0, 0, -ΔωZQ, ΔωZQ, -ΔωDQ, ΔωDQ])
    L = -1im*Ω + K

    R2eff = similar(νCPMG)
    for j = 1:length(νCPMG)
        ncyc = νCPMG[j] * T
        τ = T / 4ncyc

        Uτ = exp(L*τ)
        Uecho = Uτ * UπC * Uτ
        UechoN = Uecho^ncyc
        U = UechoN * UπH * UechoN

        relintensity = real(Mobs*U*M0)[1]
        R2eff[j] = R20 - 1/T*log(relintensity < 1e-6 ? 1e-6 : relintensity)
    end
    return CPMGData_MQ(νCPMG, R2eff, B0, T)
end

function sim_CPMG(expt::CPMGData_SQ1H, mp::modelpars_2state, ΔδH, ΔδC, R20)
    # basis for 1H SQ CPMG experiment is (Hx, Hy)
    # initial magnetisation is Hx
    # final magnetisation read-out is Hx

    # pulse program parameters
    T = expt.Trelax
    B0 = expt.B0
    νCPMG = expt.νCPMG
    kex = mp.kex
    pB = mp.pB

    ΔωH = γH * B0 * ΔδH * 1e-6

    # initial magnetisation
    M0 = [(1-pB), 0, pB, 0]
    Mobs = -[1 0 1 0]

    # set up propagators
    UπHx = [
        1 0
        0 -1
        ] |> float
    UπHy = [
        -1 0
        0 1
        ] |> float
    UπHx = kron(I(2), UπHx)
    UπHy = kron(I(2), UπHy)
    kAB = pB*kex
    kBA = (1-pB)*kex
    K = kron([
        -kAB kBA
        kAB -kBA
        ], I(2))
    Ω = [
        0 0 0 0
        0 0 0 0
        0 0 0 -ΔωH
        0 0 ΔωH 0
        ]
    L = Ω + K

    R2eff = similar(νCPMG)
    for j = 1:length(νCPMG)
        ncyc = νCPMG[j] * T
        τ = T / 4ncyc

        Uτ = exp(L*τ)
        Uecho = Uτ * UπHx * Uτ
        UechoN = Uecho^ncyc
        U = UechoN * UπHy * UechoN
        relintensity = real(Mobs*U*M0)[1]
        R2eff[j] = R20 - 1/T*log(relintensity < 1e-6 ? 1e-6 : relintensity)
    end
    return CPMGData_SQ1H(νCPMG, R2eff, B0, T)
end


function data(R2data::R2Data)
    return vec([R2data.ZQ; R2data.DQ; R2data.QQ⁻; R2data.QQ⁺])
end
function data(CPMGdata::T) where T<:CPMGData
    return CPMGdata.R2eff
end


function residuals_R2(expt::R2Data, mp, ΔδH, ΔδC, R2_ZQ, R2_DQ, R2_QQ⁻, R2_QQ⁺, B0s)
    sim = sim_R2s(mp, ΔδH, ΔδC, R2_ZQ, R2_DQ, R2_QQ⁻, R2_QQ⁺, B0s)
    return zscore.(data(expt) .- data(sim_R2s(mp, ΔδH, ΔδC, R2_ZQ, R2_DQ, R2_QQ⁻, R2_QQ⁺, B0s)))
end

function residuals_CPMG(expt::CPMGData, mp, ΔδH, ΔδC, R20)
    return zscore.(data(expt) .- data(sim_CPMG(expt, mp, ΔδH, ΔδC, R20)))
end


# """
# return residuals for particular coherence order pH
# """
# function residuals_R2(R2obs, mp, ΔδH, ΔδC, R20, pH, B0)
#     return zscore.(sim_R2(mp, ΔδH, ΔδC, R20, pH, B0) .- R2obs)
# end

# """
# return list of residuals for R2 measurements (all B0, all coherence orders)
# """
# function residuals_R2(R2obs::R2set, mp, ΔδH, ΔδC, R20s::R2set)
#     sim = sim_R2s(mp, ΔδH, ΔδC, R20s)
#     return vec([
#         zscore.(R2obs.ZQ - sim.ZQ)
#         zscore.(R2obs.DQ - sim.DQ)
#         zscore.(R2obs.QQ⁻ - sim.QQ⁻)
#         zscore.(R2obs.QQ⁺ - sim.QQ⁺)
#     ])
# end
