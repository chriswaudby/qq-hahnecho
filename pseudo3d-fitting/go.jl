using NMRTools
using Plots
using LsqFit
using Optim
using Measurements
using ColorSchemes
using LightGraphs


## input files and settings

# filenames
peaklistfilename = "data/peak.list"
spectrafilename = "data/chris_fln5_qq_160720/104/ft/test%03d.ft2"

# relaxation times (in s)
τ = [0.0011, 0.05, 0.1, 0.15]

# set size of ROI for fitting
stripwidthH = 2 * 0.025; # ppm
stripwidthC = 2 * 1.2; # ppm

# load data
spectra = loadnmr(spectrafilename)
noise = spectra[:noise]

windowC = spectra[Y, :window]
windowH = Cos²Window(0.0958464);  # not automatically parsed correctly


## includes

include("multiplet_calculations.jl")


## utility functions

findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))[2]


## load peak list

methyls = []
peakH = Vector{Float64}()
peakC = Vector{Float64}()
println("Reading peak list ($peaklistfilename):")
open(peaklistfilename) do file
    for ln in eachline(file)
        fields = split(ln)
        println("\t$(fields[1])\t$(fields[2])\t$(fields[3])")
        push!(methyls, fields[1])
        push!(peakH, parse(Float64, fields[2]))
        push!(peakC, parse(Float64, fields[3]))
    end
end
npeaks = length(methyls)
println("$npeaks peak positions loaded.")
println()


## find non-overlapping clusters
overlaps(i, j) = (abs(peakH[i]-peakH[j]) < stripwidthH) && (abs(peakC[i]-peakC[j]) < stripwidthC)
adjacencymatrix = [overlaps(i,j) for i=1:npeaks, j=1:npeaks]
clusters = connected_components(SimpleGraph(adjacencymatrix))
#clusters = [[1], [2, 3]]
nclusters = length(clusters)
println("$nclusters clusters of peaks found:")
for cluster in clusters
    println("\t$(foldr((a,b)->a*" "*b, methyls[cluster]))")
end
println()


## loop over clusters and fit
fitresults = []
for cluster in clusters

    ncluster = length(cluster)
    peaknames = methyls[cluster]
    peakgroup = foldr((a,b)->a*"-"*b, peaknames)
    peakx = peakH[cluster]
    peaky = peakC[cluster]

    # get spectrum region covering selected peaks
    xmin = minimum(peakx) - stripwidthH/2
    xmax = maximum(peakx) + stripwidthH/2
    ymin = minimum(peaky) - stripwidthC/2
    ymax = maximum(peaky) + stripwidthC/2

    # extract data ROI into arrays
    δobsH = xval(spectra[NMRTools.Between(xmin,xmax),NMRTools.Between(ymin,ymax),:])
    δobsC = yval(spectra[NMRTools.Between(xmin,xmax),NMRTools.Between(ymin,ymax),:])
    dat = data(spectra[NMRTools.Between(xmin,xmax),NMRTools.Between(ymin,ymax),:]) / noise
    nx, ny, nt = size(dat)

    # form mask for ROIs
    mask = zeros(Bool, (nx, ny))
    for i=1:ncluster
        ix = (peakx[i] - stripwidthH/2) .<= δobsH .<= (peakx[i] + stripwidthH/2)
        iy = (peaky[i] - stripwidthC/2) .<= δobsC .<= (peaky[i] + stripwidthC/2)
        mask[ix, iy] .= true
    end
    mask3 = cat([mask for i=1:nt]..., dims=3)
    maskeddat = dat[mask3]

    simmultiplet(dH, dC, R2H, R2C, S2tc, csa, J) = multiplet([dH, dC, R2H, R2C, S2tc, csa, J], δobsH, δobsC, τ, spectra[X,:bf], spectra[Y,:bf], windowH, windowC)
    # p = (ncluster x 4), p[i, :] = [R2H, R2C, S2tc, csa]
    function resid(p)
        simmultiplets = zeros(Float64, (nx, ny, nt, ncluster))
        maskedmultiplets = zeros(Float64, (length(maskeddat), ncluster))
        ysim = zeros(Float64, (nx, ny, nt))
        for i = 1:ncluster
            m = simmultiplet(peakx[i], peaky[i], p[i, 1], p[i, 2], p[i, 3], p[i, 4], 124.5)
            simmultiplets[:,:,:,i] = m
            maskedmultiplets[:, i] = m[mask3]
        end
        A = maskedmultiplets \ maskeddat
        return maskeddat - maskedmultiplets*A
    end


    # second round of fitting, including peak positions and J coupling
    function sim(p)
        simmultiplets = zeros(Float64, (nx, ny, nt, ncluster))
        maskedmultiplets = zeros(Float64, (length(maskeddat), ncluster))
        ysim = zeros(Float64, (nx, ny, nt))
        for i = 1:ncluster
            m = simmultiplet(p[i, 1], p[i, 2], p[i, 3], p[i, 4], p[i, 5], p[i, 6], p[i, 7])
            simmultiplets[:,:,:,i] = m
            maskedmultiplets[:, i] = m[mask3]
        end
        A = maskedmultiplets \ maskeddat
        for i=1:ncluster
            simmultiplets[:,:,:,i] = A[i] * simmultiplets[:,:,:,i]
        end
        return sum(simmultiplets, dims=4)[:,:,:,1]
    end
    function resid2(p)
        simmultiplets = zeros(Float64, (nx, ny, nt, ncluster))
        maskedmultiplets = zeros(Float64, (length(maskeddat), ncluster))
        ysim = zeros(Float64, (nx, ny, nt))
        for i = 1:ncluster
            m = simmultiplet(p[i, 1], p[i, 2], p[i, 3], p[i, 4], p[i, 5], p[i, 6], p[i, 7])
            simmultiplets[:,:,:,i] = m
            maskedmultiplets[:, i] = m[mask3]
        end
        A = maskedmultiplets \ maskeddat
        return maskeddat - maskedmultiplets*A
    end
    p1 = zeros(Float64, ncluster, 7)
    p1[:, 1] = peakx
    p1[:, 2] = peaky
    p1[:, 3:6] .= [50.0 40.0 10.0 30.0]
    p1[:, 7] .= 124.5

    fit = LsqFit.lmfit(p -> resid2(reshape(p, :, 7)), vec(p1), Float64[], show_trace=true, maxIter=100, x_tol=1e-3, g_tol=1e-6, autodiff=:finiteforward)

    pfit = reshape(fit.param, :, 7)
    pfiterr = reshape(stderror(fit), :, 7)
    p = pfit .± pfiterr
    ysim = sim(pfit)

    # store results
    for i=1:ncluster
        push!(fitresults, (peaknames[i], p[i,:]))
    end

    # plot results - overlay of contour plots
    clev = 10 .^LinRange(0.7, 3, 15)
    plots = []
    for i=1:nt
        # 2D overlay
        plt = heatmap(reverse(δobsH), reverse(δobsC), mask[end:-1:1,end:-1:1]')
        contour!(plt, reverse(δobsH), reverse(δobsC), dat[end:-1:1,end:-1:1,i]', c=:dodgerblue, levels=clev,
          colorbar=false, xflip=true, yflip=true, xguide="δH / ppm", yguide="δC / ppm", aspect_ratio=0.1)
        contour!(plt, reverse(δobsH), reverse(δobsC), ysim[end:-1:1,end:-1:1,i]', c=:orangered, levels=clev)
        xlims!(xmin, xmax)
        ylims!(ymin, ymax)
        push!(plots, plt)
    end
    plt = plot(plots...)
    display(plt)
    # savefig(plt, "plots/contour-$(peakgroup).pdf")

end

## print results
println()
println("methyl\tdH / ppm\tdC / ppm \tR2H / s-1\tR2C / s-1\tS2tc / ns\tcsa / ppm\tJ / Hz")
println("methyl\tdH / ppm\tdC / ppm \tR2H / s-1\tR2C / s-1\tS2tc / ns\tcsa / ppm\tJ / Hz")
for fr in fitresults
    println("$(fr[1])\t$(fr[2][1])\t$(fr[2][2])\t$(fr[2][3])\t$(fr[2][4])\t$(fr[2][5])\t$(fr[2][6])\t$(fr[2][7])")
end
println()
