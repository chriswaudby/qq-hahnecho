#!/home/nmrbox/cwaudby/julia-1.5.0-beta1/bin/julia


function proc(inputname, outputname, Δp1, Δp2)
	td = 2048
	nrelax = 14

	# input phase cycle
	ϕ1 = [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6] * 2π / 7
	ϕ2 = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2] * 2π / 3

	nphase = 21

	# sum up echo and anti-echo pathways
	ϕrx1 = -Δp1*ϕ1 - Δp2*ϕ2
	ϕrx2 =  Δp1*ϕ1 + Δp2*ϕ2

	ϕrx1 = exp.(1im * ϕrx1)
	ϕrx2 = exp.(1im * ϕrx2)

	npoints = Int(filesize(inputname)/4 - 512)
	ncomplex = Int(npoints / (td*nphase*2))
	# preallocate data (and dummy header)
	header = zeros(Float32, 512)
	y = zeros(Float32, npoints)

	# read the input file
	open(inputname) do f
	    read!(f, header)
	    read!(f, y)
	end

	y = reshape(y, td, 2, :)
	yc = y[:,1,:] + 1im * y[:,2,:]

	y = reshape(yc, td, nphase, :)
	ϕrx1 = reshape(ϕrx1, 1, nphase, 1)
	ϕrx2 = reshape(ϕrx2, 1, nphase, 1)
	y1 = sum(y .* ϕrx1, dims=2)
	y2 = sum(y .* ϕrx2, dims=2)
	y = y1 + y2
	y = reshape(y, td, ncomplex)

	# read the file
	open(outputname, "w") do f
	    write(f, header)
	    for i=1:ncomplex
			write(f, Float32.(real.(y[:,i])))
			write(f, Float32.(imag.(y[:,i])))
	    end
	end

	run(`sethdr $outputname -yN $nrelax -yT $nrelax`)
end

proc("cube.fid", "cubeQQ.fid", 3, -1)
proc("cube.fid", "cubeDQ.fid", 3, 1)
