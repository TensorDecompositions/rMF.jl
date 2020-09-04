function loadtensor(; casename = 20170911, yearstep = 1, wellsset = "01", speciesset = "20", period = 2005:yearstep:2016, datadir="tensor-data")
	global case = casename
	global wellssetid = wellsset
	global speciessetid = speciesset
	global tensoryearstep = yearstep
	global tensorperiod = period
	nt = length(period)
	data = Vector{Array{Float32,2}}(nt)
	wells = Vector{Vector{String}}(nt)
	species = Vector{Vector{String}}(nt)
	wellnameorder = Vector{Vector{String}}(nt)
	wellorder = Vector{Vector{Int64}}(nt)

	info("Load data ...")
	for (c, y) in enumerate(period)
		rMF.loaddata(casename; wellsset=wellsset, speciesset=speciesset, datemin=Date(y), datemax=Date(y+yearstep))
		wellorder[c] = rMF.wellorder
		wellnameorder[c] = convert(Array{String}, rMF.wellnameorder)
		wells[c] = convert(Array{String}, rMF.uniquewells[rMF.wellorder])
		species[c] = convert(Array{String}, rMF.uniquespecies)
		data[c] = rMF.datamatrix[rMF.wellorder, :]
		info(y)
		rMF.displayconc()
	end

	info("Print data using rMF ...")
	for y in period
		rMF.loaddata(20170911; wellsset=wellsset, speciesset=speciesset, datemin=Date(y), datemax=Date(y+yearstep), quiet=true)
		info(y)
		rMF.displayconc()
	end

	info("Print data ...")
	for (c, y) in enumerate(period)
		info(y)
		display([rMF.transposevector(["Wells"; species[c]]); wellnameorder[c] data[c][:, :]])
	end

	info("Compose tensor ...")
	tensor = Vector{Array{Float32,2}}(nt)
	wellnames = wells[end]
	for (c, y) in enumerate(period)
		tensor[c] = Array{Float32,2}(0, length(species[c]))
		for w in wellnames
			wi = findin(wells[c], [w])
			if length(wi) > 0
				tensor[c] = vcat(tensor[c], data[c][wi[1], :]')
			else
				tensor[c] = vcat(tensor[c], fill(NaN32, length(species[c]))')
			end
		end
	end

	info("Print tensor ...")
	X = Array{Float32,3}(size(tensor[end])..., length(tensor))
	for (c, y) in enumerate(period)
		info(y)
		display([rMF.transposevector(["Wells"; species[c]]); wellnames tensor[c][:, :]])
		X[:,:,c] = tensor[c]
	end

	filename = "$datadir/cr-$case-w$wellsset-s$speciesset-y$yearstep.jld"
	if !isfile(filename)
		JLD.save(filename, "X", X)
	end
	return X
end