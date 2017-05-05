"""
Perform rMF analyses
"""
function execute(range::Union{UnitRange{Int},Int}=1:maxbuckets; retries::Int=10, mixmatch::Bool=true, mixtures::Bool=true, normalize::Bool=false, scale::Bool=false, regularizationweight::Float32=convert(Float32, 0), weightinverse::Bool=false, matchwaterdeltas::Bool=false, quiet::Bool=true, clusterweights::Bool=true, convertdeltas::Bool=true, ignoreratios::Bool=false, nooutput::Bool=false)
	if sizeof(datamatrix) == 0
		warn("Execute `rMF.loaddata()` first!")
		return
	end
	concmatrix = datamatrix
	if length(ratioindex) > 0 && !ignoreratios
		concmatrix = datamatrix[:, concindex]
		ratiomatrix = datamatrix[:, ratioindex]
		nummixtures = size(concmatrix, 1)
		numconstituents = size(concmatrix, 2)
		!nooutput && info("Mixtures matrix:")
		!nooutput && display([transposevector(["Wells"; uniquespecies[concindex]]); uniquewells concmatrix])
		!nooutput && info("Ratio matrix:")
		!nooutput && display([transposevector(["Wells"; uniquespecies[ratioindex]]); uniquewells ratiomatrix])
		ratios = convert(Array{Float32, 3}, fill(NaN, nummixtures, numconstituents, numconstituents))
		for i = 1:nummixtures
			for j = 1:size(ratiocomponents, 2)
				a = ratiocomponents[1, j]
				b = ratiocomponents[2, j]
				ratios[i, a, b] = ratiomatrix[i,j] # upper triangle (not needed; added for completeness)
				ratios[i, b, a] = ratiomatrix[i,j] # lower triangle (currently used in MixMatch)
			end
		end
	else
		ratiomatrix = nothing
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
			isotopeconcentrations = MixMatch.getisotopeconcentration(deltamatrix, deltastandards, datamatrix[:, deltadependency])
			!nooutput && info("Converted delta matrix to concentrations:")
			!nooutput && display([transposevector(["Wells"; uniquespecies[deltaindex]]); uniquewells isotopeconcentrations])
			deltas = MixMatch.getisotopedelta(isotopeconcentrations, deltastandards, datamatrix[:, deltadependency])
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
	indexnan = isnan(datamatrix)
	numobservations = length(vec(datamatrix[!indexnan]))
	for numbuckets in range
		mixers[numbuckets], buckets[numbuckets], fitquality[numbuckets], robustness[numbuckets], aic[numbuckets] = NMFk.execute(concmatrix, numbuckets, retries; deltas=deltamatrix, deltaindices=deltaindices, ratios=ratiomatrix, ratioindices=ratiocomponents, mixmatch=mixmatch, normalize=normalize, scale=scale, matchwaterdeltas=matchwaterdeltas, mixtures=mixtures, quiet=quiet, regularizationweight=regularizationweight, weightinverse=weightinverse, clusterweights=clusterweights)
		mixsum = sum(mixers[numbuckets], 2)
		checkone = collect(mixsum .< 0.9) | collect(mixsum .> 1.1)
		index = find(checkone .== true)
		if length(index) > 0
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
		if casekeyword == ""
			filename = "results/$case-$numbuckets-$retries.jld"
		else
			filename = "results/$case-$casekeyword-$numbuckets-$retries.jld"
		end
		JLD.save(filename, "wells", uniquewells, "species", uniquespecies, "mixers", mixers[numbuckets], "buckets", buckets[numbuckets], "fit", fitquality[numbuckets], "robustness", robustness[numbuckets], "aic", aic[numbuckets], "regularizationweight", regularizationweight)
	end
	return
end