include("constants.jl")

"""
Parse input string s as a float, otherwise return original string.
"""
function parsestring(s)
    x = tryparse(Float64, s)
    if isnothing(x)
        x = s
    end
    return x
end

"""
Return non-comment (#) lines from input file parsed as 2D array
"""
function readfile(filename)
    lines = readlines(filename)
    dat = [parsestring.(split(strip(line))) for line in lines if line[1] ≠ '#']
    return hcat(dat...)'
end


function loadmethyls(filename)
    println("load methyls: $filename")
    return readlines(filename)
end


function loadR2(filename)
    println("load R2: $filename")

    tmp = readfile(filename)
    R2_zq = tmp[:,1:8:end] .± tmp[:,2:8:end]
    R2_dq = tmp[:,3:8:end] .± tmp[:,4:8:end]
    R2_dqq = tmp[:,5:8:end] .± tmp[:,6:8:end]
    R2_qq = tmp[:,7:8:end] .± tmp[:,8:8:end]

    return R2_zq, R2_dq, R2_dqq, R2_qq
end


function loadS2tc(filename)
    println("load S2tc: $filename")

    S2tc = readfile(filename)
    return S2tc[:,1] .± S2tc[:,2]
end


function load1Hcsa(ηHfilename, S2tc, B0_MHz)
    println("load and calculate 1H CSA:")
    println("   ηH filename: $ηHfilename")
    println("   B0 (1H, in MHz) = $B0_MHz")

    ηH = readfile(ηHfilename)
    ηH = ηH[:,1] .± ηH[:,2]

    B0 = B0_MHz * 1e6 * 2π / γH;

    ΔσH = @. 1e6*ηH/(-(4/15)*(μ₀/(4π))*ħ*γH^2*γC*P2cosβ*rCH^-3*B0*S2tc*1e-9)
    return ΔσH
end


function load13Ccsa(filename)
    println("load 13C csa: $filename")

    ΔσC = readfile(filename)
    return ΔσC[:,1] .± ΔσC[:,2]
end
