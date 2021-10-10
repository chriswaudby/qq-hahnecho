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


function loadR2(filename, B0_MHz)
    println("load R2: $filename")
    println("   B0 (1H, in MHz) = $B0_MHz")

    tmp = readfile(filename)
    B0 = B0_MHz' * 1e6 * 2π / γH;
    R2_zq = tmp[:,1:8:end] .± tmp[:,2:8:end]
    R2_dq = tmp[:,3:8:end] .± tmp[:,4:8:end]
    R2_dqq = tmp[:,5:8:end] .± tmp[:,6:8:end]
    R2_qq = tmp[:,7:8:end] .± tmp[:,8:8:end]

    #return R2set(R2_zq, R2_dq, R2_dqq, R2_qq, B0)

    data = Vector{R2Data}()
    nmethyls, ~ = size(R2_zq)
    for i=1:nmethyls
        push!(data, R2Data(R2_zq[i,:], R2_dq[i,:], R2_dqq[i,:], R2_qq[i,:], B0))
    end
    return data
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


"""
calculate and subtract CSA contributions from Hahn-echo relaxation rates (for 1H coherence order pH)
"""
function subtractCSA(R2data::Vector{R2Data}, S2tc, ΔσH, ΔσC)
    result = similar(R2data)
    for i=1:length(R2data)
        R2s = R2data[i]
        ξHcsa = 2/(3*sqrt(5)) * γH * ΔσH[i] * 1e-6 * sqrt(S2tc[i] * 1e-9)
        ξCcsa = 2/(3*sqrt(5)) * γC * ΔσC[i] * 1e-6 * sqrt(S2tc[i] * 1e-9)

        B0 = R2s.B0'
        zq = @. R2s.ZQ - (ξCcsa - ξHcsa)^2 * B0^2
        dq = @. R2s.DQ - (ξCcsa + ξHcsa)^2 * B0^2
        dqq = @. R2s.QQ⁻ - (ξCcsa - 3ξHcsa)^2 * B0^2
        qq = @. R2s.QQ⁺ - (ξCcsa + 3ξHcsa)^2 * B0^2

        result[i] = R2Data(zq, dq, dqq, qq, B0)
    end
    return result

    # ξHcsa = @. 2/(3*sqrt(5)) * γH * ΔσH * 1e-6 * sqrt(S2tc * 1e-9)
    # ξCcsa = @. 2/(3*sqrt(5)) * γC * ΔσC * 1e-6 * sqrt(S2tc * 1e-9)

    # B0 = R2s.B0
    # zq = @. R2s.ZQ - (ξCcsa - ξHcsa)^2 * B0^2
    # dq = @. R2s.DQ - (ξCcsa + ξHcsa)^2 * B0^2
    # dqq = @. R2s.QQ⁻ - (ξCcsa - 3ξHcsa)^2 * B0^2
    # qq = @. R2s.QQ⁺ - (ξCcsa + 3ξHcsa)^2 * B0^2

    # return R2set(zq, dq, dqq, qq, B0)
end


"""
load CPMG data

# Parameters:
type = :MQ or :SQ1H
Trelax = relaxation time (in seconds)
B0_MHz = 1H Larmor frequency (in MHz)
"""
function loadCPMG(filename, methyls, type, B0_MHz, Trelax; min_err=nothing)
    nmethyls = length(methyls)
    println("load CPMG data: $filename ($nmethyls methyls)")
    println("   experiment type = $(String(type))")
    println("   B₀ = $B0_MHz MHz")
    println("   Trelax = $Trelax s")
    isnothing(min_err) || println("   minimum error = $min_err s⁻¹")
    
    B0 = B0_MHz * 1e6 * 2π / γH;
    if type == :MQ
        datatype = CPMGData_MQ
    elseif type == :SQ1H
        datatype = CPMGData_SQ1H
    end
    cpmgdata = Vector{datatype}(undef, nmethyls)

    lines = readlines(filename)
    linespermethyl = length(lines) ÷ nmethyls # integer division
    νCPMG = zeros(linespermethyl-4)
    #R2 = zeros(Measurement{Float64}, linespermethyl-4, nmethyls)
    for i=1:nmethyls
        headerline = lines[1 + linespermethyl*(i-1)]
        datalines = lines[3 + linespermethyl*(i-1):linespermethyl*i-2]
        methylname = replace(headerline[2:8], 'H'=>'C') # extract assuming format Rxxx(C/H)xx, e.g. L661HD1, I743CG2, etc.
        methylindex = findfirst(methyls .== methylname)
        data = hcat([parse.(Float64, split(line)) for line in datalines]...)'
        νCPMG = data[:,1]
        err = data[:,4]
        if ~isnothing(min_err)
            @.err[err<min_err] = min_err
        end
        #R2[:,methylindex] = data[:,2] .± err
        cpmgdata[methylindex] = datatype(νCPMG, data[:,2] .± err, B0, Trelax)
    end
    #return νCPMG, R2
    return cpmgdata
end