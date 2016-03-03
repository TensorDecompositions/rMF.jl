module rMF

import NMFk
import JLD
import DataStructures

maxbuckets = 10
mixers = Array(Array{Float64, 2}, maxbuckets)
buckets = Array(Array{Float64, 2}, maxbuckets)
fitquality = Array(Float64, maxbuckets)
robustness = Array(Float64, maxbuckets)

info("rMF ...")
info("Execute `loaddata()` to get the data ...")
info("Execute `execute(5, retries=50)` to get the results for the 5 bucket case with 50 reruns ...")
info("Execute `getresults(5, retries=50)` to see the results for the 5 bucket case.")

function getresultsshort(range=1:maxbuckets; retries=10)
	for numbuckets = range
		if isfile("rmf-$case-$numbuckets-$retries.jld")
			j = JLD.load("rmf-$case-$numbuckets-$retries.jld")
			buckets[numbuckets] = j["buckets"]
			mixers[numbuckets] = j["mixers"]
			fitquality[numbuckets] = j["fit"]
			robustness[numbuckets] = j["robustness"]
			println("Buckets = $numbuckets; Best objective function = $(fitquality[numbuckets]); Robustness = $(robustness[numbuckets])")
		end
	end
end

function getresults(numbuckets=5; retries=10)
	w = false
	e = false
	try
		if !isdefined(buckets[numbuckets])
			e = true
		end
	catch
		e = true
	end
	if e 
		if isfile("rmf-$case-$numbuckets-$retries.jld")
			j = JLD.load("rmf-$case-$numbuckets-$retries.jld")
			buckets[numbuckets] = j["buckets"]
			mixers[numbuckets] = j["mixers"]
			fitquality[numbuckets] = j["fit"]
			robustness[numbuckets] = j["robustness"]
			w = true
		else
			error("Result file `rmf-$case-$numbuckets-$retries.jld` is missing ...\nExecute `execute($numbuckets)` to get the results!")
			return	
		end
	end
	info("Fit quality: $(fitquality[numbuckets])")
	info("Robustness: $(robustness[numbuckets])")
	info("Buckets:")
	display([uniquenames buckets[numbuckets]'])
	info("Match errors:")
	display([["Wells"; uniquenames]'; uniquewells concmatrix - mixers[numbuckets] * buckets[numbuckets]]')
	info("Relative match errors:")
	display([["Wells"; uniquenames]'; uniquewells (concmatrix - mixers[numbuckets] * buckets[numbuckets])./concmatrix]')
	if w
		warn("Results are loaded from external file `rmf-$case-$numbuckets-$retries.jld`...")
		warn("Execute `execute($numbuckets)` to rerun ...")
	end
end

function loaddata(timestamp="20160202")
	global case=timestamp
	dict=JLD.load("dictionary20160224.jld","dictionary")
	rawwells = readcsv("data/wells$(timestamp).csv")
	rawpz = readcsv("data/pz$(timestamp).csv")
	wells = [rawwells[2:end, 1]; rawpz[2:end, 1]]
	dates = [rawwells[2:end, 2]; rawpz[2:end, 2]]
	longnames = [rawwells[2:end, 4]; rawpz[2:end, 4]]
	names = [rawwells[2:end, 5]; rawpz[2:end, 5]]
	concs = [rawwells[2:end, 10]; rawpz[2:end, 10]]
	@assert length(longnames) == length(names)
	info("Total data record count $(length(longnames))")
	e = false
	for i in 1:length(longnames)
		if haskey(dict, longnames[i])
			names[i] = dict[longnames[i]]
		else
			e = true
			error("Variable name $(longnames[i]) is unknown!")
		end
	end
	if e
		error("Quits!")
		return
	end
	goodindices = 1:length(names)
	# remove MCOI LAOI SCI R-6i TA-53i
	goodindices = filter(i->!contains(wells[i], "MCOI"), goodindices)
	goodindices = filter(i->!contains(wells[i], "LAOI"), goodindices)
	goodindices = filter(i->!contains(wells[i], "SCI"), goodindices)
	goodindices = filter(i->!contains(wells[i], "R-6i"), goodindices)
	goodindices = filter(i->!contains(wells[i], "TA-53i"), goodindices)
	# remove ratios
	# goodindices = filter(i->!contains(longnames[i], "ratio"), goodindices)
	# goodindices = filter(i->!contains(longnames[i], "Ratio"), goodindices)
	# goodindices = filter(i->!contains(names[i], "I129"), goodindices)
	# goodindices = filter(i->!contains(names[i], "Cl36"), goodindices)
	# goodindices = filter(i->!contains(dates[i], "/2015"), goodindices) # keep only data from 2015
	wells = wells[goodindices]
	names = names[goodindices]
	dates = dates[goodindices]
	concs = concs[goodindices]
	info("Processed data record count $(length(names))")
	global uniquewells = unique(wells)
	wells2i = Dict(zip(uniquewells, 1:length(uniquewells)))
	global uniquenames = unique(names)
	name2j = Dict(zip(uniquenames, 1:length(uniquenames)))
	datacount = zeros(Int, length(uniquewells), length(uniquenames))
	concmatrix = zeros(Float64, length(uniquewells), length(uniquenames))
	for index = 1:length(wells)
		i = wells2i[wells[index]]
		j = name2j[names[index]]
		datacount[i, j] += 1
		concmatrix[i, j] += concs[index]
	end
	display([["Wells"; uniquenames]'; uniquewells concmatrix])
	display([["Wells"; uniquenames]'; uniquewells datacount])
	info("Observations per well:")
	display([uniquewells sum(datacount,2)])
	info("Observations per species:")
	display([uniquenames sum(datacount,1)'])
	global concmatrix = concmatrix ./ datacount # gives NaN if there is no data, otherwise divides by the number of results
	info("Concentration matrix:")
	display(concmatrix)
	return
end

function execute(range=1:maxbuckets; retries=10)
	if sizeof(concmatrix) == 0
		warn("Execute `loaddata()` first!")
		return
	end
	calibrationtargets = copy(concmatrix)
	for i = 1:size(calibrationtargets, 2)
		calibrationtargets[:, i] /= maximum(concmatrix[:, i]) # normalize
	end
	for numbuckets = range
		# NMFk using mixmatch 
		mixers[numbuckets], buckets[numbuckets], fitquality[numbuckets], robustness[numbuckets] = NMFk.execute(calibrationtargets, retries, numbuckets; mixmatch=true, matchdelta=false, quiet=true, regularizationweight=1e-3)
		for i = 1:size(concmatrix, 2)
			buckets[numbuckets][:, i] *= maximum(concmatrix[:, i]) # undo the normalization
		end
		println("Buckets = $numbuckets; Best objective function = $(fitquality[numbuckets]); Robustness = $(robustness[numbuckets])")
		JLD.save("rmf-$case-$numbuckets-$retries.jld", "mixers", mixers[numbuckets], "buckets", buckets[numbuckets], "fit", fitquality[numbuckets], "robustness", robustness[numbuckets])
	end
	return
end

info("Have fun ...")

end
