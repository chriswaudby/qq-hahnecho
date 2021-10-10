function plotdata(obs::T, fit=nothing) where T<:R2Data
    plt = plot(framestyle=:box,
        xguide="Magnetic field strength / T",
        yguide="Relaxation rate / s⁻¹",
        legend=:outertopright)

    scatter!(plt, obs.B0, obs.ZQ, label="ZQ")
    isnothing(fit) || plot!(plt, fit.B0, fit.ZQ, primary=false)
    scatter!(plt, obs.B0, obs.DQ, label="DQ")
    isnothing(fit) || plot!(plt, fit.B0, fit.DQ, primary=false)
    scatter!(plt, obs.B0, obs.QQ⁻, label="QQ⁻")
    isnothing(fit) || plot!(plt, fit.B0, fit.QQ⁻, primary=false)
    scatter!(plt, obs.B0, obs.QQ⁺, label="QQ⁺")
    isnothing(fit) || plot!(plt, fit.B0, fit.QQ⁺, primary=false)
    
    return plt
end


function plotdata(obs::T, fit=nothing) where T<:CPMGData
    plt = plot(framestyle=:box,
        xguide="νCPMG / Hz",
        yguide="Relaxation rate / s⁻¹",
        legend=nothing)

    scatter!(plt, obs.νCPMG, obs.R2eff)
    isnothing(fit) || plot!(plt, fit.νCPMG, fit.R2eff, primary=false)
    
    return plt
end


function plotdata(obs::T, fits=nothing) where T<:Vector{CPMGData}
    plt = plot(framestyle=:box,
        xguide="νCPMG / Hz",
        yguide="Relaxation rate / s⁻¹",
        legend=:outertopright)

    if isnothing(fits)
        fits = fill(nothing,length(obs))
    end
    for (ob, fit) in zip(obs, fits)
        if isa(ob, CPMGData_MQ)
            type = "MQ"
        elseif isa(ob, CPMGData_SQ1H)
            type = "1H SQ"
        else
            type = "unknown CPMG type"
        end
        label = "$(round(ob.B0,digits=1)) T, $type"
        scatter!(plt, ob.νCPMG, ob.R2eff, label=label)
        isnothing(fit) || plot!(plt, fit.νCPMG, fit.R2eff, primary=false)
    end
    
    return plt
end
