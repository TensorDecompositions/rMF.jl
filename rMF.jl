module rMF

import MixMatch
import NMFk
import JLD
import DataStructures
import SpatialAnalysis
import PyCall
import Gadfly
import Colors
import PyPlot

@PyCall.pyimport matplotlib.patheffects as PathEffects

maxbuckets = 10
mixers = Array(Array{Float64, 2}, maxbuckets)
buckets = Array(Array{Float64, 2}, maxbuckets)
fitquality = Array(Float64, maxbuckets)
robustness = Array(Float64, maxbuckets)
deltadependancy = Int[]
dict = Dict()

#=
dict = DataStructures.OrderedDict{AbstractString,AbstractString}(
	"Chromium"=>"Cr",
	"Chromium-53/52"=>"δ53Cr",
	"Bromide"=>"Br-",
	"Chloride"=>"Cl-",
	"Chlorate"=>"ClO3",
	"Perchlorate"=>"ClO4",
	"Chlorine-36/Chlorine Ratio (e-15)"=>"Cl36",
	"Tritium"=>"3H",
	"Deuterium Ratio"=>"δ2H",
	"Oxygen-18/Oxygen-16 Ratio"=>"δ18O",
	"Nitrate-Nitrite as Nitrogen"=>"NO3",
	"Nitrogen-15/Nitrogen-14 Ratio(NO3)"=>"δ15N",
	"Oxygen-18/Oxygen-16 Ratio from Nitrate"=>"δ18O-NO3",
	"Sulfate"=>"SO4",
	"Sulfur-34/Sulfur-32 Ratio (SO4)"=>"δ34S-SO4",
	"Oxygen-18/Oxygen-16 Ratio from SO4"=>"δ18O-SO4",
	"Iodine-129/Iodine Ratio (e-15)"=>"I129",
	"Fraction Modern Carbon (de-normalized)"=>"f14C",
	"Dioxane[1,4-]"=>"Dioxane",
	"Acetaminophen"=>"Acetam",
	"Caffeine"=>"Caffe",
	"Sulfamethoxazole"=>"Sulfame")
JLD.save("dictionary20160202ordered.jld","dictionary",dict)
=#

info("rMF (Robust Matrix Factorization)")
info("")
info("Use `rMF.loaddata()` to get data")
info("Use `rMF.execute(5, retries=50)` to get the results for the 5 bucket case with 50 reruns")
info("Use `rMF.getresults(5, retries=50)` to see the results for the 5 bucket case.")
info("Use `rMF.getresultsshort(5:8, retries=50)` to see the results for bucket cases 5 to 8.")

function getresultsshort(range=1:maxbuckets; retries=10)
	for numbuckets = range
		if isfile("results/$case-$numbuckets-$retries.jld")
			j = JLD.load("results/$case-$numbuckets-$retries.jld")
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
		if isfile("results/$case-$numbuckets-$retries.jld")
			j = JLD.load("results/$case-$numbuckets-$retries.jld")
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
				if isfile("dictionary$(case).jld")
					dictold=JLD.load("dictionary$(case).jld","dictionary")
					remap = DataStructures.OrderedDict(zip(collect(keys(dictold)), 1:length(uniquespecies)))
				else
					remap = DataStructures.OrderedDict(zip(uniquespecies_long, 1:length(uniquespecies)))
				end
			end
			for i in keys(order)
				order[i] = remap[i]
			end
			buckets[numbuckets] = buckets[numbuckets][:,collect(values(order))]
			w = true
		else
			error("Result file `results/$case-$numbuckets-$retries.jld` is missing ...\nExecute `rMF.execute($numbuckets)` to get the results!")
			return	
		end
	end
	info("Fit quality: $(fitquality[numbuckets])")
	info("Robustness: $(robustness[numbuckets])")

	wells2i = Dict(zip(uniquewells, 1:length(uniquewells)))
	if isfile("orderedwells-WE.dat")
		wellnameorder = readdlm("orderedwells-WE.dat")
		wellorder = zeros(Int, length(wellnameorder))
		for i in 1:length(wellorder)
			if haskey(wells2i, wellnameorder[i])
				wellorder[i] = wells2i[wellnameorder[i]]
			end
		end
		wellmissing = wellorder .== 0
		indexmissing = find(wellmissing)
		wellavailable = wellorder .!= 0
		indexavailale = find(wellavailable)
		if length(indexmissing) > 0 && length(indexavailale) != 0
			warn("Well(s) $(wellnameorder[indexmissing]) is not accounted!")
		else
			wellorder = 1:length(uniquewells)
			wellnameorder = uniquewells
		end
	else
		wellorder = 1:length(uniquewells)
		wellnameorder = uniquewells		
	end

	if length(deltaindex) > 0
		numconstituents = size(datamatrix, 2)
		numdeltas = length(deltaindex)
		numconc = numconstituents - numdeltas
		H_conc = buckets[numbuckets][:,1:numconc]
		H_deltas = buckets[numbuckets][:,numconc+1:end]
		predictions = similar(datamatrix)
		predictions[:, concindex] = mixers[numbuckets] * H_conc
		predictions[:, deltaindex] = MixMatch.computedeltas(mixers[numbuckets], H_conc, H_deltas, deltadependancy)
	else
		predictions = mixers[numbuckets] * buckets[numbuckets]
	end
	if length(ratioindex) > 0
		info("Compute ratios:")
		p = similar(datamatrix)
		p[:, concindex] = predictions
		for i = 1:size(predictions, 1)
			for j = 1:size(ratiocomponents, 2)
				a = ratiocomponents[1, j]
				b = ratiocomponents[2, j]
				p[i, ratioindex[j]] = p[i, a] / p[i, b]
			end
		end
		predictions = p
	end

	errors = datamatrix - predictions
	relerrors = errors ./ datamatrix

	f = open("results/$case-data.dat", "w")
	writedlm(f, [["Wells"; uniquespecies]'; wellnameorder datamatrix[wellorder, :]]')
	close(f)

	info("Predictions:")
	display([["Wells"; uniquespecies]'; wellnameorder predictions[wellorder, :]]')
	f = open("results/$case-$numbuckets-$retries-predictions.dat", "w")
	writedlm(f, [["Wells"; uniquespecies]'; wellnameorder predictions[wellorder, :]]')
	close(f)
	info("Match errors:")
	display([["Wells"; uniquespecies]'; wellnameorder errors[wellorder, :]]')
	f = open("results/$case-$numbuckets-$retries-errors.dat", "w")
	writedlm(f, [["Wells"; uniquespecies]'; wellnameorder errors[wellorder, :]]')
	close(f)
	info("Relative match errors:")
	display([["Wells"; uniquespecies]'; wellnameorder relerrors[wellorder, :]]')
	f = open("results/$case-$numbuckets-$retries-relerrors.dat", "w")
	writedlm(f, [["Wells"; uniquespecies]'; wellnameorder relerrors[wellorder, :]]')
	close(f)
	info("Max/min Species in Buckets:")
	maxs = maximum(buckets[numbuckets], 1)
	mins = minimum(buckets[numbuckets], 1)
	display([uniquespecies[dataindex] maxs' mins'])
	info("Mixers:")
	display([uniquewells mixers[numbuckets]])
	info("Buckets:")
	display([uniquespecies[dataindex] buckets[numbuckets]'])
	info("Normalized buckets:")
	for i=1:length(maxs)
		if maxs[i] == mins[i]
			mins[i] = 0
		end
	end
	nbuckets = (buckets[numbuckets] .- mins) ./ (maxs - mins)
	display([uniquespecies[dataindex] nbuckets'])
	source_weight = sum(nbuckets, 2)
	source_index = sortperm(collect(source_weight), rev=true)
	info("Sorted normalized buckets:")
	sbuckets = nbuckets[source_index,:]
	display([uniquespecies[dataindex] sbuckets'])
	info("Sorted buckets:")
	s2buckets = buckets[numbuckets][source_index,:]
	display([uniquespecies[dataindex] s2buckets'])
	# sbuckets[sbuckets.<1e-6] = 1e-6
	gbucket = Gadfly.spy(sbuckets', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
				Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"),
				Gadfly.Theme(default_point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
				Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red"))))
	# filename, format = Mads.setimagefileformat(filename, format)
	filename = "results/$case-$numbuckets-$retries-buckets.svg"
	Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
	filename = "results/$case-$numbuckets-$retries-buckets.png"
	Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)

	info("Sorted mixers:")
	smixers = mixers[numbuckets][wellorder, source_index]
	display([wellnameorder smixers])
	gmixers = Gadfly.spy(smixers, Gadfly.Scale.y_discrete(labels = i->wellnameorder[i]), Gadfly.Scale.x_discrete,
				Gadfly.Guide.YLabel("Wells"), Gadfly.Guide.XLabel("Sources"),
				Gadfly.Theme(default_point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
				Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red"))))
	# filename, format = Mads.setimagefileformat(filename, format)
	filename = "results/$case-$numbuckets-$retries-mixers.svg"
	Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)
	filename = "results/$case-$numbuckets-$retries-mixers.png"
	Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)

	return
# remove deep screens
	goodindices = 1:length(uniquewells)
	goodindices = filter(i->!contains(uniquewells[i], "_2"), goodindices)
	wn = uniquewells[goodindices,:]
	wn = map(i->replace(wn[i], "_1", ""), 1:length(wn))
	wc = wellcoord[goodindices, :]
	if length(wc) == 0
		warn("No well coordinates")
		return
	end
	grid_size = 100
	mm_grid_size = grid_size * 5
	minx = floor(minimum(wc, 1)[1] / mm_grid_size - .5) * mm_grid_size
	maxx = ceil(maximum(wc, 1)[1] / mm_grid_size + .5) * mm_grid_size
	miny = floor(minimum(wc, 1)[2] / mm_grid_size - .5) * mm_grid_size
	maxy = ceil(maximum(wc, 1)[2] / mm_grid_size + .5) * mm_grid_size
	xi = minx:grid_size:maxx
	yi = miny:grid_size:maxy
	zi = Array(Float64, length(xi), length(yi))
	wm = mixers[numbuckets][goodindices,source_index]
	for s = 1:numbuckets
		z = wm[:, s] 
		z[z.<1e-6] = 1e-6
		z = log10(z)
		@assert length(wc[:,1]) == length(z)
		@assert length(wc[:,1]) == length(unique(wc[:,1])) || length(wc[:,2]) == length(unique(wc[:,2]))
		zi = SpatialAnalysis.linear_interpolation(wc[:,1], wc[:,2], z, xi, yi)
		zmin = minimum(zi)
		zmin = -3
		zmax = maximum(zi)
		dz = ( zmax - zmin ) / 2
		plt = PyPlot.figure(figsize=(9, 4))
		PyPlot.imshow(zi, origin="lower", extent=[minx, maxx, miny, maxy], vmin=zmin-dz, vmax=zmax+dz, cmap="jet")
		PyPlot.colorbar(shrink=0.5, cmap="jet")
		PyPlot.clim(zmin, zmax)
		PyPlot.title("Source $s")
		PyPlot.scatter(wc[:,1], wc[:,2], marker="o", color="black", c=z, s=70, cmap="jet")
		PyPlot.scatter(wc[:,1], wc[:,2], marker="o", color="white", c=z, s=68, cmap="jet")
		PyPlot.clim(zmin, zmax)
		for i = 1:length(wn)
			PyPlot.annotate(wn[i], xy=(wc[i,1], wc[i,2]), xytext=(-2, 2), fontsize=8, color="white", textcoords="offset points", ha="right", va="bottom", path_effects=[PathEffects.withStroke(linewidth=1, foreground="black")])
		end
		PyPlot.yticks([538500, 539500], ["538500", "539500"], rotation="vertical")
		PyPlot.tight_layout()
		PyPlot.savefig("results/$case-$numbuckets-$retries-source-$s.png")
		PyPlot.close()
	end
	selectedspecies =[1,3,5,7]
	swm = map(i->rMF.mixers[5][goodindices,i] * rMF.buckets[5][i,selectedspecies], source_index)
	for b = 1:numbuckets
		for s = 1:length(selectedspecies)
			z = swm[b][:, s]
			zmin = minimum(z)
			zmax = maximum(z)
			current_species = uniquespecies[selectedspecies[s]]
			info("Source $b: $current_species")
			info("Min = $(zmin)")
			info("Max = $(zmax)")
			# z[z.<1e-6] = 1e-6
			z = log10(z)
			@assert length(wc[:,1]) == length(z)
			@assert length(wc[:,1]) == length(unique(wc[:,1])) || length(wc[:,2]) == length(unique(wc[:,2]))
			zi = SpatialAnalysis.linear_interpolation(wc[:,1], wc[:,2], z, xi, yi)
			zmin = minimum(zi)
			zmin = -3
			zmax = maximum(zi)
			if zmin > zmax
				zmin = minimum(zi)
			end
			dz = ( zmax - zmin ) / 2
			plt = PyPlot.figure(figsize=(9, 4))
			PyPlot.imshow(zi, origin="lower", extent=[minx, maxx, miny, maxy], vmin=zmin-dz, vmax=zmax+dz, cmap="jet")
			PyPlot.colorbar(shrink=0.5, cmap="jet")
			PyPlot.clim(zmin, zmax)
			PyPlot.title("Source $b: $current_species")
			PyPlot.scatter(wc[:,1], wc[:,2], marker="o", color="black", c=z, s=70, cmap="jet")
			PyPlot.scatter(wc[:,1], wc[:,2], marker="o", color="white", c=z, s=68, cmap="jet")
			PyPlot.clim(zmin, zmax)
			for i = 1:length(wn)
				PyPlot.annotate(wn[i], xy=(wc[i,1], wc[i,2]), xytext=(-2, 2), fontsize=8, color="white", textcoords="offset points", ha="right", va="bottom", path_effects=[PathEffects.withStroke(linewidth=1, foreground="black")])
			end
			PyPlot.yticks([538500, 539500], ["538500", "539500"], rotation="vertical")
			PyPlot.tight_layout()
			PyPlot.savefig("results/$case-$numbuckets-$retries-source-$b-$current_species.png")
			PyPlot.close()
		end
	end
	if w
		warn("Results are loaded from external file `results/$case-$numbuckets-$retries.jld`...")
		warn("Execute `rMF.execute($numbuckets)` to rerun ...")
	end
end

function loaddata(casename::AbstractString)
	global case = casename
	global ratioindex = Int[]
	global deltas = nothing
	global deltaindex = Int[]
	if casename == "test23delta"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "δA"]
		global uniquespecies_long = uniquespecies
		global datamatrix = [[1., 0.1] [0.1, 1.] [0.1, 1.]]
		global deltaindex = Int[3]
		global deltadependancy = Int[1]
		global dataindex = collect(1:size(datamatrix,2))
		global concindex = setdiff(dataindex, deltaindex)
		global wellcoord = [[0., 0.] [0., 100.]]
		info("Concentration matrix:")
		display([["Wells"; uniquespecies]'; uniquewells datamatrix])
		return
	end
	if casename == "test23ratio"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "A/B"]
		global uniquespecies_long = uniquespecies
		global datamatrix = [[1., NaN] [NaN, 2.] [1., 1.]]
		global ratioindex = Int[3]
		global ratiocomponents = Int[1, 2]
		global concindex = setdiff(collect(1:size(datamatrix,2)), ratioindex)
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.]]
		info("Concentration matrix:")
		display([["Wells"; uniquespecies]'; uniquewells datamatrix])
		return
	end
	if casename == "test23"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "C"]
		global uniquespecies_long = uniquespecies
		global datamatrix = [[1., 1.] [2., 4.] [2., 2.]]
		global concindex = collect(1:size(datamatrix,2))
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.]]
		info("Concentration matrix:")
		display([["Wells"; uniquespecies]'; uniquewells datamatrix])
		return
	end
	rawdata = readcsv("data/" * casename * ".csv")
	global uniquewells = rawdata[2:end, 1]
	global uniquespecies = rawdata[1, 2:end]'
	global uniquespecies_long = uniquespecies
	global datamatrix = rawdata[2:end, 2:end]
	global concindex = collect(1:size(datamatrix,2))
	global dataindex = concindex
	info("Species ($(length(uniquespecies)))")
	display(uniquespecies)
	info("Wells ($(length(uniquewells)))")
	display(uniquewells)
	info("Concentration matrix:")
	datamatrix[datamatrix .== " "] = NaN
	datamatrix[datamatrix .== ""] = NaN
	datamatrix = Array{Float32}(datamatrix)
	display([["Wells"; uniquespecies]'; uniquewells datamatrix])
	return
end

function loaddata(probstamp::Int64=20160102, probkey::AbstractString="")
	dictstamp = probstamp
	global case = probstamp
	global dict = JLD.load("dictionary$(dictstamp)ordered.jld","dictionary")
	# mixes = JLD.load("mixtures$(timestamp).jld","mixtures")
	rawwells = readcsv("data/wells$(probstamp).csv")[2:end,[1,2,4,5,10]]
	rawpz = readcsv("data/pz$(probstamp).csv")[2:end,[1,2,4,5,10]]
	rawdata = [rawwells; rawpz]
	wells = rawdata[:, 1]
	dates = rawdata[:, 2]
	longnames = rawdata[:, 3]
	names = rawdata[:, 4]
	concs = rawdata[:, 5]
	@assert length(longnames) == length(names)
	@assert length(names) == length(concs)
	info("Total data record count $(length(longnames))")
	# set the "cool" (short) variable names
	for i in 1:length(longnames)
		if haskey(dict, longnames[i])
			names[i] = dict[longnames[i]]
		end
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
	longnames = longnames[goodindices]
	dates = dates[goodindices]
	concs = concs[goodindices]
	info("Processed data record count $(length(names))")
	global uniquewells = unique(wells)
	wells2i = Dict(zip(uniquewells, 1:length(uniquewells)))
	dnames = collect(values(dict))
	global uniquespecies = dnames
	global uniquespecies_long = collect(keys(dict))
	name2j = Dict(zip(uniquespecies, 1:length(uniquespecies)))
	datacount = zeros(Int, length(uniquewells), length(uniquespecies))
	datamatrix = zeros(Float64, length(uniquewells), length(uniquespecies))
	for index = 1:length(wells)
		i = wells2i[wells[index]]
		if haskey(name2j, names[index])
			j = name2j[names[index]]
			datacount[i, j] += 1
			datamatrix[i, j] += concs[index]
		end
	end
	info("Species ($(length(uniquespecies)))")
	display(uniquespecies)
	info("Wells ($(length(uniquewells)))")
	display(uniquewells)
	info("Observations count:")
	display([["Wells"; uniquespecies]'; uniquewells datacount])
	info("Observations per well:")
	display([uniquewells sum(datacount,2)])
	info("Observations per species:")
	display([uniquespecies sum(datacount,1)'])
	global datamatrix = datamatrix ./ datacount # gives NaN if there is no data, otherwise divides by the number of results
	info("Concentration matrix:")
	display([["Wells"; uniquespecies]'; uniquewells datamatrix])
	global concindex = collect(1:size(datamatrix,2))
	coord, coordheader = readdlm("data/coord$(probstamp).dat", header=true)
	global wellcoord = Array(Float64, length(uniquewells), 2)
	for index = 1:length(uniquewells)
		i = indexin([uniquewells[index]], coord[:,1])[1]
		if i == 0
			warn("Coordinates for well $(uniquewells[index]) are missing!")
		else
			wellcoord[index,1] = coord[i,2]
			wellcoord[index,2] = coord[i,3]
		end
	end
	info("Check species in the dictionary ...")
	uniquelongnames = unique(longnames)
	for i in 1:length(uniquelongnames)
		if !haskey(dict, uniquelongnames[i])
			warn("Species name $(uniquelongnames[i]) in the data set is not defined in the dictionary!")
		end
	end
	info("Check species in the data set ...")
	for i in 1:length(uniquespecies_long)
		if indexin([uniquespecies_long[i]], uniquelongnames)[1] == 0
			warn("Species name $(uniquespecies_long[i]) defined in the dictionary is missing in the data set!")
		end
	end
	return
end

"""
Perform rMF analyses
"""
function execute(range=1:maxbuckets; retries::Int=10, mixmatch::Bool=true, normalize::Bool=false, mixtures::Bool=true, regularizationweight::Number=0, matchwaterdeltas::Bool=false, quiet::Bool=true)
	if sizeof(datamatrix) == 0
		warn("Execute `rMF.loaddata()` first!")
		return
	end
	if length(ratioindex) > 0
		nummixtures = size(datamatrix, 1)
		numconstituents = size(datamatrix, 2)
		concmatrix = datamatrix[:, concindex]
		ratiomatrix = datamatrix[:, ratioindex]
		info("Mixtures matrix:")
		display([["Wells"; uniquespecies[concindex]]'; uniquewells concmatrix])
		info("Ratio matrix:")
		display([["Wells"; uniquespecies[ratioindex]]'; uniquewells ratiomatrix])
		ratios = fill(NaN, nummixtures, numconstituents, numconstituents)
		for i = 1:nummixtures
			for j = 1:size(ratiocomponents, 2)
				a = ratiocomponents[1, j]
				b = ratiocomponents[2, j]
				ratios[i, a, b] = ratiomatrix[i]
			end
		end
	else
		concmatrix = datamatrix
		ratios = nothing
	end
	if length(deltaindex) > 0
		nummixtures = size(datamatrix, 1)
		numconstituents = size(datamatrix, 2)
		concmatrix = datamatrix[:, concindex]
		deltamatrix = datamatrix[:, deltaindex]
		info("Mixtures matrix:")
		display([["Wells"; uniquespecies[concindex]]'; uniquewells concmatrix])
		info("Delta matrix:")
		display([["Wells"; uniquespecies[deltaindex]]'; uniquewells deltamatrix])
	else
		concmatrix = datamatrix
		deltamatrix = nothing
	end
	for numbuckets = range
		# NMFk using mixmatch 
		mixers[numbuckets], buckets[numbuckets], fitquality[numbuckets], robustness[numbuckets] = NMFk.execute(concmatrix, retries, numbuckets; deltas=deltamatrix, deltaindices=deltadependancy, ratios=ratios, mixmatch=mixmatch, normalize=normalize, matchwaterdeltas=matchwaterdeltas, mixtures=mixtures, quiet=quiet, regularizationweight=regularizationweight)
		mixsum = sum(mixers[numbuckets], 2)
		index = find(mixsum .> 1.1) | find(mixsum .< 0.9)
		if length(index) > 0
			warn("The mixers do not add to 1")
			@show sum(mixers[numbuckets], 2)
			display(mixers[numbuckets])
		end
		println("Buckets = $numbuckets; Best objective function = $(fitquality[numbuckets]); Robustness = $(robustness[numbuckets])")
		JLD.save("results/$case-$numbuckets-$retries.jld", "wells", uniquewells, "species", uniquespecies, "mixers", mixers[numbuckets], "buckets", buckets[numbuckets], "fit", fitquality[numbuckets], "robustness", robustness[numbuckets])
	end
	return
end

info("")
info("Use `rMF.loaddata()` to load different data sets:")
info("rMF.loaddata(20151202) - original fingerprint data set")
info("rMF.loaddata(20160102) - original fingerprint data set without pz wells")
info("rMF.loaddata(20160202) - new fingerprint data set with pz wells")
info("rMF.loaddata(20160302) - new fingerprint data set without pz wells")
info("""rMF.loaddata("rdx-20160721")""")
info("")
info("Use `rMF.execute()` to perform rMF analyses:")
info("rMF.execute(2:4, retries=100)")
info("")
info("Have fun ...")

# cd(Pkg.dir("rMF") * "/AquiferMixing")

# loaddata()

end
