module rMF

import NMFk
import JLD
import DataStructures
using Gadfly
using Colors

maxbuckets = 10
mixers = Array(Array{Float64, 2}, maxbuckets)
buckets = Array(Array{Float64, 2}, maxbuckets)
fitquality = Array(Float64, maxbuckets)
robustness = Array(Float64, maxbuckets)

#=
dict = DataStructures.OrderedDict{AbstractString,AbstractString}(
	"Chromium"=>"Cr",
	"Chromium-53/52 Ratio"=>"δ53Cr",
	"Chloride"=>"Cl-",
	"Chlorate"=>"ClO3",
	"Perchlorate"=>"ClO4",
	"Chlorine-36"=>"Cl36",
	"Tritium"=>"3H",
	"Deuterium Ratio"=>"δ2H",
	"Oxygen-18/Oxygen-16 Ratio"=>"δ18O",
	"Nitrate-Nitrite as Nitrogen"=>"NO3",
	"Nitrogen-15/Nitrogen-14 Ratio(NO3)"=>"δ15N",
	"Oxygen-18/Oxygen-16 Ratio from Nitrate"=>"δ18O-NO3",
	"Sulfate"=>"SO4",
	"Sulfur-34/Sulfur-32 Ratio (SO4)"=>"δ34S-SO4",
	"Oxygen-18/Oxygen-16 Ratio from SO4"=>"δ18O-SO4",
	"Iodine-129"=>"I129",
	"Dioxane[1,4-]"=>"Dioxane",
	"Acetaminophen"=>"Acetam",
	"Caffeine"=>"Caffe",
	"Sulfamethoxazole"=>"Sulfame")
JLD.save("dictionary20160224ordered.jld","dictionary",dict)
=#

info("rMF ...")
info("Execute `rMF.loaddata()` to get the data ...")
info("Execute `rMF.execute(5, retries=50)` to get the results for the 5 bucket case with 50 reruns ...")
info("Execute `rMF.getresults(5, retries=50)` to see the results for the 5 bucket case.")
info("Execute `rMF.getresultsshort(5:8, retries=50)` to see the results for bucket cases 5 to 8.")

function getresultsshort(range=1:maxbuckets; retries=10)
	for numbuckets = range
		if isfile("results/rmf-$case-$numbuckets-$retries.jld")
			j = JLD.load("results/rmf-$case-$numbuckets-$retries.jld")
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
		if isfile("results/rmf-$case-$numbuckets-$retries.jld")
			j = JLD.load("results/rmf-$case-$numbuckets-$retries.jld")
			buckets[numbuckets] = j["buckets"]
			mixers[numbuckets] = j["mixers"]
			fitquality[numbuckets] = j["fit"]
			robustness[numbuckets] = j["robustness"]
			order = DataStructures.OrderedDict(zip(uniquespecies_long, 1:length(uniquespecies)))
			if haskey(j, "uniquewells") && haskey(j, "uniquespecies")
				wells = j["uniquewells"]
				species = j["uniquespecies"]
				remap = DataStructures.OrderedDict(zip(species, 1:length(uniquespecies)))
			else
				dictold=JLD.load("dictionary$(case).jld","dictionary")
				remap = DataStructures.OrderedDict(zip(collect(keys(dictold)), 1:length(uniquespecies)))
			end
			for i in keys(order)
				order[i] = remap[i]
			end
			buckets[numbuckets] = buckets[numbuckets][:,collect(values(order))]
			w = true
		else
			error("Result file `results/rmf-$case-$numbuckets-$retries.jld` is missing ...\nExecute `execute($numbuckets)` to get the results!")
			return	
		end
	end
	info("Fit quality: $(fitquality[numbuckets])")
	info("Robustness: $(robustness[numbuckets])")
	info("Match errors:")
	display([["Wells"; uniquespecies]'; uniquewells concmatrix - mixers[numbuckets] * buckets[numbuckets]]')
	info("Relative match errors:")
	display([["Wells"; uniquespecies]'; uniquewells (concmatrix - mixers[numbuckets] * buckets[numbuckets])./concmatrix]')
	info("Max/min Species:")
	maxs = maximum(buckets[numbuckets], 1)
	mins = minimum(buckets[numbuckets], 1)
	display([uniquespecies maxs' mins'])
	info("Buckets:")
	display([uniquespecies buckets[numbuckets]'])
	info("Normalized buckets:")
	nbuckets = ( buckets[numbuckets] .- mins ) ./ (maxs - mins)
	display([uniquespecies nbuckets'])
	s = sum(nbuckets, 2)
	i = sortperm(collect(s), rev=true)
	info("Sorted normalized buckets:")
	sbuckets = nbuckets[i,:]
	display([uniquespecies sbuckets'])
	info("Sorted buckets:")
	s2buckets = buckets[numbuckets][i,:]
	display([uniquespecies s2buckets'])
	# sbuckets[sbuckets.<1e-6] = 1e-6
	gbucket = Gadfly.spy(sbuckets', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
				Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"),
				Gadfly.Theme(default_point_size=20pt, major_label_font_size=14pt, minor_label_font_size=12pt, key_title_font_size=16pt, key_label_font_size=12pt),
				Gadfly.Scale.ContinuousColorScale(Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red"))))
	# filename, format = Mads.setimagefileformat(filename, format)
	filename = "results/rmf-$case-$numbuckets-$retries-buckets.svg"
	Gadfly.draw(Gadfly.SVG(filename,6inch,12inch), gbucket)
	filename = "results/rmf-$case-$numbuckets-$retries-buckets.png"
	Gadfly.draw(Gadfly.PNG(filename,6inch,12inch), gbucket)
	if w
		warn("Results are loaded from external file `results/rmf-$case-$numbuckets-$retries.jld`...")
		warn("Execute `execute($numbuckets)` to rerun ...")
	end
end

function loaddata(timestamp="20160102")
	global case=timestamp
	dict=JLD.load("dictionary$(timestamp)ordered.jld","dictionary")
	rawwells = readcsv("data/wells$(timestamp).csv")[2:end,[1,2,4,5,10]]
	rawpz = readcsv("data/pz$(timestamp).csv")[2:end,[1,2,4,5,10]]
	rawdata = [rawwells; rawpz]
	wells = rawdata[:, 1]
	dates = rawdata[:, 2]
	longnames = rawdata[:, 3]
	names = rawdata[:, 4]
	concs = rawdata[:, 5]
	@assert length(longnames) == length(names)
	@assert length(names) == length(concs)
	info("Total data record count $(length(longnames))")
	# check the variable names
	e = false
	for i in 1:length(longnames)
		if haskey(dict, longnames[i])
			names[i] = dict[longnames[i]]
		else
			e = true
			warn("Variable name $(longnames[i]) is unknown!")
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
	dnames = collect(values(dict))
	@assert length(dnames) == length(unique(names))
	global uniquespecies = dnames
	global uniquespecies_long = collect(keys(dict))
	name2j = Dict(zip(uniquespecies, 1:length(uniquespecies)))
	datacount = zeros(Int, length(uniquewells), length(uniquespecies))
	concmatrix = zeros(Float64, length(uniquewells), length(uniquespecies))
	for index = 1:length(wells)
		i = wells2i[wells[index]]
		j = name2j[names[index]]
		datacount[i, j] += 1
		concmatrix[i, j] += concs[index]
	end
	info("Observations count:")
	display([["Wells"; uniquespecies]'; uniquewells datacount])
	info("Observations per well:")
	display([uniquewells sum(datacount,2)])
	info("Observations per species:")
	display([uniquespecies sum(datacount,1)'])
	global concmatrix = concmatrix ./ datacount # gives NaN if there is no data, otherwise divides by the number of results
	info("Concentration matrix:")
	display([["Wells"; uniquespecies]'; uniquewells concmatrix])
	coord, coordheader = readdlm("data/coord$(timestamp).dat", header=true)
	for index = 1:length(uniquewells)
		if indexin([uniquewells[index]], coord[:,1])[1] == 0
			warn("Coordinates for well $(uniquewells[index]) are missing!")
		end
	end
	global wellcoord = coord
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
		JLD.save("results/rmf-$case-$numbuckets-$retries.jld", "wells", uniquewells, "species", uniquespecies, "mixers", mixers[numbuckets], "buckets", buckets[numbuckets], "fit", fitquality[numbuckets], "robustness", robustness[numbuckets])
	end
	return
end

info("Have fun ...")
info("""loaddata(timestamp="20151202") - original fingerprint data set""")
info("""loaddata(timestamp="20160102") - original fingerprint data set without pz wells""")
info("""loaddata(timestamp="20160202") - new fingerprint data set without pz wells""")

loaddata()

end
