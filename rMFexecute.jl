"Test rMF functions"
function test()
	include(joinpath(rmfdir, "test", "runtests.jl"))
end

"""
Perform rMF analyses
"""
function execute(range::Union{Range{Int},Int}=1:maxbuckets; retries::Int=10, normalize::Bool=false, scale::Bool=false, regularizationweight::Float32=convert(Float32, 0), weightinverse::Bool=false, quiet::Bool=true, clusterweights::Bool=false, convertdeltas::Bool=true, ignoreratios::Bool=false, nooutput::Bool=false, mixture::Symbol=:mixmatch, method::Symbol=:ipopt, save::Bool=true, resultdir::AbstractString="results", tol=1e-3, kw...)
	if !isdefined(rMF, :datamatrix) || sizeof(datamatrix) == 0
		warn("rMF problem is not defined; execute `rMF.loaddata()` first!")
		return
	end
	if !isdir(resultdir)
		mkdir(resultdir)
	end
	if casekeyword == ""
		casefilename = joinpath(resultdir, "$case")
	else
		casefilename = joinpath(resultdir, "$case-$casekeyword")
	end
	concmatrix = datamatrix
	numberofratios = length(ratioindex)
	if numberofratios > 0 && !ignoreratios
		if length(size(ratiocomponents)) == 2
			nr, nc = size(ratiocomponents)
		else
			nc = length(ratiocomponents)
			nr = 1
		end
		if numberofratios != nr && nc != 2 && length(ratiocomponents[ratiocomponents.==0]) != 0
			warn("Ratio data is corrupted; ratios will be ignored!")
			global ratiomatrix = Array{Float32}(0, 0)
			global ratiocomponents = Array{Int}(0, 0)
		else
			concmatrix = datamatrix[:, concindex]
			ratiomatrix = datamatrix[:, ratioindex]
			nummixtures = size(concmatrix, 1)
			numconstituents = size(concmatrix, 2)
			!nooutput && info("Mixtures matrix:")
			!nooutput && display([transposevector(["Wells"; uniquespecies[concindex]]); uniquewells concmatrix])
			!nooutput && info("Ratio matrix:")
			!nooutput && display([transposevector(["Wells"; uniquespecies[ratioindex]]); uniquewells ratiomatrix])
		end
	else
		global ratiomatrix = Array{Float32}(0, 0)
		global ratiocomponents = Array{Int}(0, 0)
	end
	if length(deltaindex) > 0
		concmatrix = datamatrix[:, concindex]
		deltamatrix = datamatrix[:, deltaindex]
		deltaindices = indexin(deltadependency, concindex)
		!nooutput && info("Mixtures matrix:")
		!nooutput && display([transposevector(["Wells"; uniquespecies[concindex]]); uniquewells concmatrix])
		!nooutput && info("Delta matrix:")
		!nooutput && display([transposevector(["Wells"; uniquespecies[deltaindex]]); uniquewells deltamatrix])
		if convertdeltas
			isotopeconcentrations = NMFk.getisotopeconcentration(deltamatrix, deltastandards, datamatrix[:, deltadependency])
			!nooutput && info("Converted delta matrix to concentrations:")
			!nooutput && display([transposevector(["Wells"; uniquespecies[deltaindex]]); uniquewells isotopeconcentrations])
			deltas = NMFk.getisotopedelta(isotopeconcentrations, deltastandards, datamatrix[:, deltadependency])
			!nooutput && info("Converted stable isotope concentrations back to deltas (test):")
			!nooutput && display([transposevector(["Wells"; uniquespecies[deltaindex]]); uniquewells deltas])
			concmatrix = copy(datamatrix)
			concmatrix[:, deltaindex] = isotopeconcentrations
			!nooutput && info("Concentration matrix:")
			!nooutput && display([transposevector(["Wells"; uniquespecies]); uniquewells concmatrix])
			deltaindices = deltadependency
			deltamatrix = Array{Float32}(0, 0)
		end
	else
		deltaindices = deltadependency
		deltamatrix = Array{Float32}(0, 0)
	end
	indexnan = isnan.(datamatrix)
	numobservations = length(vec(datamatrix[.!indexnan]))
	for numbuckets in range
		mixers[numbuckets], buckets[numbuckets], fitquality[numbuckets], robustness[numbuckets], aic[numbuckets] = NMFk.execute(concmatrix, numbuckets, retries; deltas=deltamatrix, deltaindices=deltaindices, ratios=ratiomatrix, ratioindices=ratiocomponents, normalize=normalize, scale=scale, quiet=quiet, regularizationweight=regularizationweight, weightinverse=weightinverse, clusterweights=clusterweights, mixture=mixture, method=method, tol=tol, casefilename=casefilename, kw...)
		mixsum = sum(mixers[numbuckets], 2)
		checkone = collect(mixsum .< 0.9) .| collect(mixsum .> 1.1)
		index = find(checkone .== true)
		if length(index) > 0 && mixture != :null
			warn("Mixer matrix rows do not add to 1")
			display(mixsum)
			@show mixers[numbuckets]
			@show buckets[numbuckets]
		end
		if convertdeltas && isdefined(:deltastandards)
			deltas = MixMatch.getisotopedelta(buckets[numbuckets][:, deltaindex], deltastandards, buckets[numbuckets][:, deltadependency])
			buckets[numbuckets][:, deltaindex] = deltas
		end
		if length(range) == 1
			info("Estimated buckets:")
			display([uniquespecies[dataindex] buckets[numbuckets]'])
			if sizeof(truebucket) > 0
				info("True buckets:")
				if sizeof(truedeltas) > 0
					display([uniquespecies[dataindex] hcat(truebucket, truedeltas)'])
				else
					display([uniquespecies[dataindex] truebucket'])
				end
			end
			info("Estimated mixers:")
			display([uniquewells mixers[numbuckets]])
			if sizeof(truemixer) > 0
				info("True mixers:")
				display([uniquewells truemixer])
			end
		end
		if save
			filename = "$casefilename-$numbuckets-$retries.jld"
			JLD.save(filename, "wells", convert(Array{String,1}, uniquewells), "species", convert(Array{String,1}, uniquespecies), "mixers", mixers[numbuckets], "buckets", buckets[numbuckets], "fit", fitquality[numbuckets], "robustness", robustness[numbuckets], "aic", aic[numbuckets], "regularizationweight", regularizationweight)
		end
	end
	return
end