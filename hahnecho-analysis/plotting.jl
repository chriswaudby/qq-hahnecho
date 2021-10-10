function plotχ2(hval, cval, χ2, βZQ, βDQ, βDQQ, βQQ)
    χ2min = minimum(vec(χ2))
    ξHmin = hval[argmin(χ2)[1]]
    ξCmin = cval[argmin(χ2)[2]]

    clev = log10.(χ2min .+ [2.3, 6.7])
    plt = plot(framestyle=:box, legend=nothing, xguide="ξH", yguide="ξC", title="$(methyls[resi])", xlims=(-.3, .3), ylims=(-.3, .3))

    plot!(plt, hval,   hval .+ Measurements.value(sqrt(βZQ-(ξCcsa[resi]-ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βZQ-(ξCcsa[resi]-ξHcsa[resi])^2)), c=1)
    plot!(plt, hval,   hval .- Measurements.value(sqrt(βZQ-(ξCcsa[resi]-ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βZQ-(ξCcsa[resi]-ξHcsa[resi])^2)), c=1)
    plot!(plt, hval,  -hval .+ Measurements.value(sqrt(βDQ-(ξCcsa[resi]+ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βDQ-(ξCcsa[resi]+ξHcsa[resi])^2)), c=2)
    plot!(plt, hval,  -hval .- Measurements.value(sqrt(βDQ-(ξCcsa[resi]+ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βDQ-(ξCcsa[resi]+ξHcsa[resi])^2)), c=2)
    plot!(plt, hval,  3hval .+ Measurements.value(sqrt(βDQQ-(ξCcsa[resi]-3ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βDQQ-(ξCcsa[resi]-3ξHcsa[resi])^2)), c=3)
    plot!(plt, hval,  3hval .- Measurements.value(sqrt(βDQQ-(ξCcsa[resi]-3ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βDQQ-(ξCcsa[resi]-3ξHcsa[resi])^2)), c=3)
    plot!(plt, hval, -3hval .+ Measurements.value(sqrt(βQQ-(ξCcsa[resi]+3ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βQQ-(ξCcsa[resi]+3ξHcsa[resi])^2)), c=4)
    plot!(plt, hval, -3hval .- Measurements.value(sqrt(βQQ-(ξCcsa[resi]+3ξHcsa[resi])^2)),
        ribbon=0hval .+ Measurements.uncertainty(sqrt(βQQ-(ξCcsa[resi]+3ξHcsa[resi])^2)), c=4)

    contour!(plt, hval, cval, log10.(χ2'), levels=clev, c=:red)

    return plt
end

function plotlinefits(ξH, ξC, resi)
    βZQ = Measurements.value((ξCcsa[resi]-ξHcsa[resi])^2 + (ξC-ξH)^2)
    βDQ = Measurements.value((ξCcsa[resi]+ξHcsa[resi])^2 + (ξC+ξH)^2)
    βDQQ = Measurements.value((ξCcsa[resi]-3ξHcsa[resi])^2 + (ξC-3ξH)^2)
    βQQ = Measurements.value((ξCcsa[resi]+3ξHcsa[resi])^2 + (ξC+3ξH)^2)
    
    plt = plot(frame=:box, legend=:topleft, xguide="B₀² / T²", yguide="Relaxation rate / s⁻¹")
    c0 = [10.0]
    
    scatter!(plt, B0², zq[resi,:],label="ZQ")
    fit = LsqFit.lmfit(c -> resid(B0², zq[resi,:], βZQ, c), c0, Float64[])
    c = fit.param[1]
    B2 = [0; B0²]
    plot!(plt, B2, line(B2, βZQ, c),primary=false)

    scatter!(plt, B0², dq[resi,:],label="DQ")
    fit = LsqFit.lmfit(c -> resid(B0², dq[resi,:], βDQ, c), c0, Float64[])
    c = fit.param[1]
    B2 = [0; B0²]
    plot!(plt, B2, line(B2, βDQ, c),primary=false)

    scatter!(plt, B0², dqq[resi,:],label="DQQ")
    fit = LsqFit.lmfit(c -> resid(B0², dqq[resi,:], βDQQ, c), c0, Float64[])
    c = fit.param[1]
    B2 = [0; B0²]
    plot!(plt, B2, line(B2, βDQQ, c),primary=false)

    scatter!(plt, B0², qq[resi,:],label="QQ")
    fit = LsqFit.lmfit(c -> resid(B0², qq[resi,:], βQQ, c), c0, Float64[])
    c = fit.param[1]
    B2 = [0; B0²]
    plot!(plt, B2, line(B2, βQQ, c),primary=false)
    ylims!(0, ylims()[2])
    return plt
end