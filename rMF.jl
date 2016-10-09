module rMF

import MixMatch
import NMFk
import JLD
import StatsBase
import DataStructures
import SpatialAnalysis
import PyCall
import Gadfly
import Compose
import Colors
import PyPlot
import Images
import Compat
import Compat.AbstractString
include(ENV["HOME"] * "/Julia/madsdisplay.jl")

@PyCall.pyimport matplotlib.patheffects as PathEffects

function transposevector(a)
	reshape(a, 1, length(a))
end

function transposematrix(a)
	permutedims(a, (2, 1))
end

maxbuckets = 10
case = ""
casekeyword = ""
mixers = Array(Array{Float64, 2}, maxbuckets)
buckets = Array(Array{Float64, 2}, maxbuckets)
fitquality = Array(Float64, maxbuckets)
robustness = Array(Float64, maxbuckets)
deltadependency = Int[]
dict_species = Dict()
uniquewells = []
uniquespecies = []

#=
dict_species = DataStructures.OrderedDict{AbstractString,AbstractString}(
	"Chromium"=>"Cr",
	"Chromium-53/52"=>"δ53Cr",
	"Bromide"=>"Br-",
	"Chloride"=>"Cl-",
	"Chlorate"=>"ClO3",
	"Perchlorate"=>"ClO4",
	"Chlorine-36/Chlorine Ratio (e-15)"=>"r36Cl",
	"Tritium"=>"3H",
	"Deuterium Ratio"=>"δ2H",
	"Oxygen-18/Oxygen-16 Ratio"=>"δ18O",
	"Nitrate-Nitrite as Nitrogen"=>"NO3",
	"Nitrogen-15/Nitrogen-14 Ratio(NO3)"=>"δ15N",
	"Oxygen-18/Oxygen-16 Ratio from Nitrate"=>"δ18O-NO3",
	"Sulfate"=>"SO4",
	"Sulfur-34/Sulfur-32 Ratio (SO4)"=>"δ34S-SO4",
	"Oxygen-18/Oxygen-16 Ratio from SO4"=>"δ18O-SO4",
	"Iodine-129/Iodine Ratio (e-15)"=>"rI129",
	"Fraction Modern Carbon (de-normalized)"=>"f14C",
	"Dioxane[1,4-]"=>"Dioxane",
	"Acetaminophen"=>"Acetam",
	"Caffeine"=>"Caffe",
	"Sulfamethoxazole"=>"Sulfame")
JLD.save("data/cr-species.jld", "species", dict_species)

mixtures = DataStructures.OrderedDict{Any,Any}(
	"Cr"=>Any[],
	"δ53Cr"=>ASCIIString["Cr"],
	"Cl-"=>Any[],
	"ClO3"=>Any[],
	"ClO4"=>Any[],
	"δ36Cl"=>ASCIIString["Cl-"],
	"3H"=>Any[],
	"δ2H"=>Any[],
	"δ18O"=>Any[],
	"NO3"=>Any[],
	"δ15N"=>ASCIIString["NO3"],
	"δ18O-NO3"=>ASCIIString["NO3"],
	"SO4"=>Any[],
	"δ34S-SO4"=>ASCIIString["SO4"],
	"δ18O-SO4"=>ASCIIString["SO4"],
	"δI129"=>Any[],
	"Dioxane"=>Any[],
	"Acetam"=>Any[],
	"Caffe"=>Any[],
	"Sulfame"=>Any[])
JLD.save("data/cr-stable-isotope-mixtures.jld", "mixtures", mixtures)
=#

info("rMF (Robust Matrix Factorization)")
info("")
info("Use `rMF.loaddata()` to get data")
info("Use `rMF.execute(5, retries=50)` to compute the results for the 5 bucket case with 50 reruns")
info("Use `rMF.getresults(5, retries=50)` to get the results for the 5 bucket case.")
info("Use `rMF.getresults(5:6, retries=50; brief=true)` to see the results for bucket cases 5 to 8.")

function getresults(range::Union{UnitRange{Int},Int}=1:maxbuckets, keyword::AbstractString=""; retries=10, brief::Bool=false)
	if keyword != ""
		if case != "" && !contains(keyword, case)
			casestring = case * "-" * keyword
		else
			casestring = keyword
		end
	else
		if case != "" 
			if casekeyword != ""
				casestring = case * "-" * casekeyword
			else
				casestring = case
			end
		else
			warn("Problem is not defined; use rMF.loaddata() first")
			return
		end
	end
	for numbuckets = range
		filename = "results/$(casestring)-$numbuckets-$retries.jld"
		if isfile(filename)
			j = JLD.load(filename)
			buckets[numbuckets] = j["buckets"]
			mixers[numbuckets] = j["mixers"]
			fitquality[numbuckets] = j["fit"]
			robustness[numbuckets] = j["robustness"]
			order = DataStructures.OrderedDict(zip(uniquespecies_long, 1:length(uniquespecies_long)))
			if haskey(j, "uniquewells") && haskey(j, "uniquespecies")
				wells = j["uniquewells"]
				species = j["uniquespecies"]
				remap = DataStructures.OrderedDict(zip(species, 1:length(species)))
			else
				if isfile("data/cr-species-order-$(case).jld")
					dictold=JLD.load("data/cr-species-order-$(case).jld","dictionary")
					remap = DataStructures.OrderedDict(zip(collect(keys(dictold)), 1:length(uniquespecies)))
				else
					remap = DataStructures.OrderedDict(zip(uniquespecies_long, 1:length(uniquespecies_long)))
				end
			end
			for i in keys(order)
				order[i] = remap[i]
			end
			if length(ratioindex) > 0
				index_list = setdiff(collect(values(order)), ratioindex)
			else				
				index_list = collect(values(order))
			end
			buckets[numbuckets] = buckets[numbuckets][:, index_list]
		else
			error("Result file `$(filename)` is missing ...\nExecute `rMF.execute($numbuckets)` to get the results!")
			continue
		end
		
		wellorder, wellnameorder = getwellorder()

		numwells = size(datamatrix, 1)
		numconstituents = size(datamatrix, 2)
		@assert numwells == length(wellnameorder)
		@assert numconstituents == length(uniquespecies)

		orderedbuckets = similar(buckets[numbuckets])
		global spredictions = Array(Array{Float64, 2}, numbuckets)
		if length(deltaindex) > 0
			numdeltas = length(deltaindex)
			numconc = numconstituents - numdeltas
			H_conc = buckets[numbuckets][:,1:numconc]
			H_deltas = buckets[numbuckets][:,numconc+1:end]
			orderedbuckets[:, concindex] = H_conc
			orderedbuckets[:, deltaindex] = H_deltas
			global predictions = similar(datamatrix)
			predictions[:, concindex] = mixers[numbuckets] * H_conc
			predictions[:, deltaindex] = MixMatch.computedeltas(mixers[numbuckets], H_conc, H_deltas, indexin(deltadependency, concindex))
			tpredictions = zeros(size(datamatrix))
			for i = 1:numbuckets
				spredictions[i] = similar(datamatrix)
				spredictions[i][:, concindex] = mixers[numbuckets][:,i:i] * H_conc[i:i,:]
				spredictions[i][:, deltaindex] = MixMatch.computedeltas(reshape(mixers[numbuckets][:,i:i], numwells, 1), H_conc[i:i,:], H_deltas[i:i,:], indexin(deltadependency, concindex), compute_contributions=true)
				tpredictions = tpredictions + spredictions[i]
			end
			tpredictions[:, deltaindex] .*= predictions[:, deltaindex] ./ tpredictions[:, deltaindex]
			@Base.Test.test_approx_eq_eps maximum(abs(predictions .- tpredictions)) 0 1e-3
		else
			predictions = mixers[numbuckets] * buckets[numbuckets]
			orderedbuckets = buckets[numbuckets]
			for i = 1:numbuckets
				spredictions[i] = similar(datamatrix)
				spredictions[i] = mixers[numbuckets][:,i:i] * buckets[numbuckets][i:i,:]
			end
		end
		if length(ratioindex) > 0
			info("Compute ratios:")
			p = similar(datamatrix)
			p[:, concindex] = predictions
			for i = 1:numwells
				for j = 1:size(ratiocomponents, 2)
					a = ratiocomponents[1, j]
					b = ratiocomponents[2, j]
					p[i, ratioindex[j]] = p[i, a] / p[i, b]
				end
			end
			predictions = p
			for i = 1:numbuckets
				p = similar(datamatrix)
				p[:, concindex] = spredictions[i]
				for i = 1:numwells
					for j = 1:size(ratiocomponents, 2)
						a = ratiocomponents[1, j]
						b = ratiocomponents[2, j]
						p[i, ratioindex[j]] = p[i, a] / p[i, b]
					end
				end
				spredictions[i] = p
			end
		end

		indexnan = isnan(datamatrix)
		d = copy(datamatrix)
		d[indexnan] = 0
		errors = d - predictions
		errors[indexnan] = 0.
		of = sum(errors.^2)
		errors[indexnan] = NaN
		relerrors = errors ./ d
		relerrors[indexnan] = NaN

		if brief
			println("Buckets = $numbuckets; Best objective function = $(fitquality[numbuckets]) (check = $(of); regularization penalty = $(sum(log(1.+abs(orderedbuckets)).^2)/numbuckets)); Robustness = $(robustness[numbuckets])")
			continue
		else
			info("Fit quality: $(fitquality[numbuckets]) (check = $(of)) (regularization penalty = $(sum(log(1.+abs(orderedbuckets)).^2)/numbuckets))")
			if of - fitquality[numbuckets] > 1e-1
				warn("Objective function test fails!")
			end
			info("Robustness: $(robustness[numbuckets])")
		end

		f = open("results/$(casestring)-data.dat", "w")
		writedlm(f, transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder datamatrix[wellorder, :]]))
		close(f)

		source_weight = sum(mixers[numbuckets], 1)
		source_index = sortperm(vec(collect(source_weight)), rev=true)

		info("Predictions:")
		display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder predictions[wellorder, :]]))
		f = open("results/$(casestring)-$numbuckets-$retries-predictions.dat", "w")
		writedlm(f, transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder predictions[wellorder, :]]))
		close(f)
		
		info("Predictions for each bucket:")
		for i = 1:numbuckets
			info("Predictions for bucket #$(i)")
			display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder spredictions[source_index[i]][wellorder, :]]))
			f = open("results/$(casestring)-$numbuckets-$retries-bucket-$i-predictions.dat", "w")
			writedlm(f, transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder spredictions[source_index[i]][wellorder, :]]))
		end

		MArows = 5
		MAcols = Int(ceil(numwells/5))
		MA = Array(Compose.Context, (MArows, MAcols))
		i = 1
		for w in wellorder
			b = abs(hcat(map(i->collect(spredictions[i][w,:]), 1:numbuckets)...)) ./ abs(predictions[w:w, :]')
			b = b ./ maximum(b, 2)
			MA[i] = Gadfly.render(Gadfly.spy(b[:,source_index], Gadfly.Guide.title(wellnameorder[i]),
					Gadfly.Guide.xticks(orientation=:horizontal), Gadfly.Guide.yticks(orientation=:horizontal),
					Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete(),
					Gadfly.Guide.YLabel("", orientation=:vertical), Gadfly.Guide.XLabel("", orientation=:horizontal), 
					Gadfly.Guide.colorkey(""),
					Gadfly.Theme(key_position=:none, major_label_font_size=10Gadfly.pt, minor_label_font_size=8Gadfly.pt, key_title_font_size=10Gadfly.pt, key_label_font_size=8Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1)))
			i += 1
		end
		for w in numwells+1:MArows*MAcols
			MA[w] = Compose.context()
		end
		gs = Gadfly.gridstack(MA)
		filename = "results/$(casestring)-$numbuckets-$retries-wellmixtures.svg"
		Gadfly.draw(Gadfly.SVG(filename, MArows * (1.6Gadfly.inch + numbuckets * 0.1Gadfly.inch), MAcols * (1.6Gadfly.inch + numconstituents * 0.1Gadfly.inch)), gs)
		filename = "results/$(casestring)-$numbuckets-$retries-wellmixtures.png"
		Gadfly.draw(Gadfly.PNG(filename, MArows * (1.6Gadfly.inch + numbuckets * 0.1Gadfly.inch), MAcols * (1.6Gadfly.inch + numconstituents * 0.1Gadfly.inch)), gs)

		info("Match errors:")
		display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder errors[wellorder, :]]))
		f = open("results/$(casestring)-$numbuckets-$retries-errors.dat", "w")
		writedlm(f, transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder errors[wellorder, :]]))
		close(f)
		
		info("Match error statistics:")
		verrrors = vec(errors[!isnan(errors)])
		numobservations = length(verrrors)
		dof = numobservations - numbuckets
		sml = dof + numobservations * (log(fitquality[numbuckets]/dof) / 2 + 1.837877)
		aic = sml + 2 * numbuckets
		println("AIC: $(aic)")
		println("Mean: $(mean(verrrors))")
		println("Variance: $(var(verrrors))")
		println("Standard Deviation: $(std(verrrors))")
		println("Maximum: $(maximum(verrrors))")
		println("Minimum: $(minimum(verrrors))")
		println("Skewness: $(StatsBase.skewness(verrrors))")
		println("Kurtosis: $(StatsBase.kurtosis(verrrors))")
		f = open("results/$(casestring)-$numbuckets-$retries-stats.dat", "w")
		println(f, "Number of buckets: $(numbuckets)")
		println(f, "* Degrees of freedom: $(dof)")
		println(f, "* Fit quality: $(fitquality[numbuckets])")
		println(f, "* Robustness: $(robustness[numbuckets])")
		println(f, "* AIC: $(aic)")
		println(f, "* Error stats: $(aic)")
		println(f, "  - Mean: $(mean(verrrors))")
		println(f, "  - Variance: $(var(verrrors))")
		println(f, "  - Standard Deviation: $(std(verrrors))")
		println(f, "  - Maximum: $(maximum(verrrors))")
		println(f, "  - Minimum: $(minimum(verrrors))")
		println(f, "  - Skewness: $(StatsBase.skewness(verrrors))")
		println(f, "  - Kurtosis: $(StatsBase.kurtosis(verrrors))")
		close(f)

		info("Histogram of the estimation errors:")
		g = Gadfly.plot(x=verrrors, Gadfly.Geom.histogram())
		filename = "results/$(casestring)-$numbuckets-$retries-error_histogram.png"
		Gadfly.draw(Gadfly.PNG(filename, 6Gadfly.inch, 4Gadfly.inch), g)
		madsdisplay(filename)

		indmaxerror = ind2sub(size(errors), indmax(abs(errors)))
		info("The largest absolute match error is for $(wellnameorder[indmaxerror[1]]) / $(uniquespecies[indmaxerror[2]]).")
		println("Observation: $(datamatrix[indmaxerror...])")
		println("Prediction: $(predictions[indmaxerror...])")
		println("Error: $(errors[indmaxerror...])")
		println("Relative error: $(relerrors[indmaxerror...])")
		display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder[indmaxerror[1]] errors[indmaxerror[1]:indmaxerror[1], :]]))
		
		info("Relative match errors:")
		display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder relerrors[wellorder, :]]))
		f = open("results/$(casestring)-$numbuckets-$retries-relerrors.dat", "w")
		writedlm(f, transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder relerrors[wellorder, :]]))
		close(f)
		
		indmaxerror = ind2sub(size(relerrors), indmax(abs(relerrors)))
		info("The largest absolute relative match error is for $(wellnameorder[indmaxerror[1]]) / $(uniquespecies[indmaxerror[2]]).")
		println("Observation: $(datamatrix[indmaxerror...])")
		println("Prediction: $(predictions[indmaxerror...])")
		println("Error: $(errors[indmaxerror...])")
		println("Relative error: $(relerrors[indmaxerror...])")
		display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder[indmaxerror[1]] relerrors[indmaxerror[1]:indmaxerror[1], :]]))
		
		info("Max/min Species in Buckets per species:")
		maxs1 = maximum(orderedbuckets, 1)
		mins1 = minimum(orderedbuckets, 1)
		display([uniquespecies[dataindex] maxs1' mins1'])
		
		info("Max/min Species in Buckets per buckets:")
		maxs2 = maximum(orderedbuckets, 2)
		mins2 = minimum(orderedbuckets, 2)
		display([maxs2 mins2])
		
		info("Mixers:")
		display([uniquewells mixers[numbuckets]])
		
		info("Buckets:")
		display([uniquespecies[dataindex] orderedbuckets[source_index,:]'])
		f = open("results/$(casestring)-$numbuckets-$retries-buckets.dat", "w")
		writedlm(f, [uniquespecies[dataindex] orderedbuckets[source_index,:]'])
		close(f)

		info("Normalized buckets:")
		for i=1:length(maxs1)
			if maxs1[i] == mins1[i]
				mins1[i] = 0
			end
		end
		for i=1:length(maxs2)
			if maxs2[i] == mins2[i]
				mins2[i] = 0
			end
		end
		nbuckets = (orderedbuckets .- mins1) ./ (maxs1 - mins1)
		display([uniquespecies[dataindex] nbuckets[source_index,:]'])
		
		info("Ordered buckets normalized by (max/min species) the overall species dominance:")
		s1buckets = nbuckets[source_index,:]
		display([uniquespecies[dataindex] s1buckets'])
		
		info("Ordered buckets normalized by (max/min buckets) species dominance within each bucket:")
		n2buckets = (orderedbuckets .- mins2) ./ (maxs2 - mins2)
		s2buckets = n2buckets[source_index,:]
		display([uniquespecies[dataindex] s2buckets'])
		bucketimpact = Array(Float64, numbuckets, numconstituents)
		for s in 1:numconstituents
			for i = 1:numbuckets
				# bucketimpact[i, s] = sum(abs((spredictions[i][:, s])./predictions[:,s]))
				bucketimpact[i, s] = sum(abs((spredictions[i][:, s])))
			end
		end
		bucketimpactwells = Array(Float64, numbuckets, numwells)
		for w in 1:numwells
			for i = 1:numbuckets
				# bucketimpactwells[i, w] = sum(abs((spredictions[i][w, :])./predictions[w,:]))
				bucketimpactwells[i, w] = sum(abs((spredictions[i][w, :])))
			end
		end
		
		info("Ordered buckets to capture the overall impact on the species concentrations:")
		display([uniquespecies[dataindex] bucketimpact[source_index, dataindex]'])
		
		info("Max/min Species in model predictions for each bucket:")
		maxm = maximum(bucketimpact, 1)
		minm = minimum(bucketimpact, 1)
		display([uniquespecies[dataindex] maxm' minm'])
		for i=1:length(maxm)
			if maxm[i] == minm[i]
				minm[i] = 0
			end
		end
		maxiw = maximum(bucketimpactwells, 1)
		miniw = minimum(bucketimpactwells, 1)
		display([wellnameorder maxiw' miniw'])
		for i=1:length(maxiw)
			if maxiw[i] == miniw[i]
				miniw[i] = 0
			end
		end
		bucketimpactwells[source_index, wellorder] = (bucketimpactwells .- miniw) ./ (maxiw - miniw)
		gmixers = Gadfly.spy(bucketimpactwells', Gadfly.Scale.y_discrete(labels = i->wellnameorder[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Wells"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(default_point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpactwells.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpactwells.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)

		info("Ordered buckets normalized to capture the overall impact on the species concentrations:")
		bucketimpact[source_index, dataindex] = (bucketimpact .- minm) ./ (maxm - minm)
		display([uniquespecies[dataindex] bucketimpact'])

		gbucket = Gadfly.spy(bucketimpact', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(default_point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpact.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpact.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)

		# s1buckets[s1buckets.<1e-6] = 1e-6
		gbucket = Gadfly.spy(s1buckets', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(default_point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets1.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets1.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)

		gbucket = Gadfly.spy(s2buckets', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(default_point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets2.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets2.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)

		info("Ordered buckets:")
		display([uniquespecies[dataindex] orderedbuckets[source_index,:]'])

		info("Ordered mixers:")
		smixers = mixers[numbuckets][wellorder, source_index]
		smixers[smixers .< 0] = 0
		smixers[smixers .> 1] = 1
		display([wellnameorder smixers])

		gmixers = Gadfly.spy(smixers, Gadfly.Scale.y_discrete(labels = i->wellnameorder[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Wells"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(default_point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-mixers.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)
		filename = "results/$(casestring)-$numbuckets-$retries-mixers.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)

		if !isdefined(:wellcoord)
			continue
		end

		# remove deep screens
		goodindices = 1:length(uniquewells)
		goodindices = filter(i->!contains(uniquewells[i], "_2"), goodindices)
		wn = uniquewells[goodindices,:]
		wn = map(i->replace(wn[i], "_1", ""), 1:length(wn))
		wc = wellcoord[goodindices, :]
		if length(wc) == 0
			warn("No well coordinates")
			continue
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
			dz = (zmax - zmin) / 2
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
			PyPlot.savefig("results/$(casestring)-$numbuckets-$retries-source-$s.png")
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
				dz = (zmax - zmin) / 2
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
				PyPlot.savefig("results/$(casestring)-$numbuckets-$retries-source-$b-$current_species.png")
				PyPlot.close()
			end
		end
	end
end

function loaddata(casename::AbstractString, keyword::AbstractString="")
	global case = casename
	global casekeyword = keyword
	global mixers = Array(Array{Float64, 2}, maxbuckets)
	global buckets = Array(Array{Float64, 2}, maxbuckets)
	global fitquality = Array(Float64, maxbuckets)
	global robustness = Array(Float64, maxbuckets)
	global ratioindex = Int[]
	global deltaindex = Int[]
	global deltadependency = Int[]
	if casename == "test23delta"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["δA", "A", "B"]
		global uniquespecies_long = uniquespecies
		global datamatrix = convert(Array{Float32, 2}, [[0.1, 1.] [1., 0.1] [0.1, 1.]])
		global deltaindex = Int[1]
		global deltadependency = Int[2]
		global dataindex = collect(1:size(datamatrix, 2))
		global concindex = setdiff(dataindex, deltaindex)
		global wellcoord = [[0., 0.] [0., 100.]]
		info("Concentration matrix:")
		display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test23delta2"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "δA"]
		global uniquespecies_long = uniquespecies
		global datamatrix = convert(Array{Float32, 2}, [[1., 0.1] [0.1, 1.] [0.1, 1.]])
		global deltaindex = Int[3]
		global deltadependency = Int[1]
		global dataindex = collect(1:size(datamatrix, 2))
		global concindex = setdiff(dataindex, deltaindex)
		global wellcoord = [[0., 0.] [0., 100.]]
		info("Concentration matrix:")
		display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test23ratio"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "A/B"]
		global uniquespecies_long = uniquespecies
		global datamatrix = convert(Array{Float32, 2}, [[NaN, 2] [1, NaN] [1., 2.]])
		global ratioindex = Int[3]
		global ratiocomponents = Int[1, 2]
		global concindex = setdiff(collect(1:size(datamatrix,2)), ratioindex)
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.]]
		info("Concentration matrix:")
		display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test23"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "C"]
		global uniquespecies_long = uniquespecies
		global datamatrix = convert(Array{Float32,2}, [[1., 1.] [2., 4.] [2., 2.]])
		global concindex = collect(1:size(datamatrix,2))
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.]]
		info("Concentration matrix:")
		display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	rawdata = readcsv("data/" * casename * ".csv")
	rawdata[rawdata .== " "] = NaN
	rawdata[rawdata .== ""] = NaN
	global uniquewells = rawdata[2:end, 1]
	global uniquespecies = rawdata[1, 2:end]'
	global uniquespecies_long = uniquespecies
	global datamatrix = convert(Array{Float32,2}, rawdata[2:end, 2:end])
	global concindex = collect(1:size(datamatrix,2))
	global dataindex = concindex
	info("Species ($(length(uniquespecies)))")
	display(uniquespecies)
	info("Wells ($(length(uniquewells)))")
	display(uniquewells)
	info("Concentration matrix:")
	display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
	return
end

function getwellorder()
	wells2i = Dict(zip(uniquewells, 1:length(uniquewells)))
	if isfile("data/cr-well-order-WE.dat")
		wellnameorder = readdlm("data/cr-well-order-WE.dat")
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
		wellorder = wellorder[indexavailale]
		wellnameorder = wellnameorder[indexavailale]
	else
		warn("data/cr-well-order-WE.dat is missing!")
		wellorder = 1:length(uniquewells)
		wellnameorder = uniquewells
	end
	if length(wellorder) == 0
		wellorder = 1:length(uniquewells)
		wellnameorder = uniquewells
	end
	return wellorder, wellnameorder
end

function displayconc(name::AbstractString)
	wellorder, wellnameorder = getwellorder()
	if name == ""
		display([transposevector(["Wells"; uniquespecies]); wellnameorder datamatrix[wellorder,:]])
	else
		i = findin(uniquespecies, [name])
		if length(i) > 0
			j = i[1]
			display([transposevector(["Wells"; uniquespecies[j]]); wellnameorder datamatrix[wellorder,j]])
		else
			i = findin(wellnameorder, [name])
			if length(i) > 0
				j = i[1]
				display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder[j] datamatrix[wellorder[j]:wellorder[j],:]]))
			end
		end
	end
end

function loaddata(probstamp::Int64=20160102, keyword::AbstractString=""; wellsset::AbstractString="", speciesset::AbstractString="")
	casestring = keyword
	if wellsset != ""
		if casestring != ""
			casestring = casestring * "-w" * wellsset
		else
			casestring = "w" * wellsset
		end
	end
	if speciesset != ""
		if casestring != ""
			casestring = casestring * "-s" * speciesset
		else
			casestring = "s" * speciesset
		end
	end
	global casekeyword = casestring
	global case = "cr-$probstamp"
	global mixers = Array(Array{Float64, 2}, maxbuckets)
	global buckets = Array(Array{Float64, 2}, maxbuckets)
	global fitquality = Array(Float64, maxbuckets)
	global robustness = Array(Float64, maxbuckets)
	global ratioindex = Int[]
	filename = "data/cr-species.jld"
	if isfile(filename)
		global dict_species = JLD.load(filename, "species")
		global uniquespecies = collect(values(dict_species))
		global uniquespecies_long = collect(keys(dict_species))
		if speciesset != ""
			filename = "data/cr-species-set$(speciesset).txt"
			if isfile(filename)
				ss = readdlm(filename)
				sd = setdiff(ss, uniquespecies)
				if length(sd) > 0
					error("There are species in $filename missing!")
					display(sd)
					return
				end
				sd = setdiff(uniquespecies, ss)
				if length(sd) > 0
					warn("The following species will be removed!")
					display(sd)
					ind = findin(uniquespecies, ss)
					global uniquespecies = uniquespecies[ind]
					global uniquespecies_long = uniquespecies_long[ind]
				end
				dnames = ss
			else
				error("$filename is missing!")
				return
			end
		end
	else
		error("$filename is missing!")
		return
	end
	filename = "data/cr-stable-isotope-mixtures.jld"
	if isfile(filename)
		isotope_mixtures = JLD.load(filename, "mixtures")
		deltaindex = Int[]
		deltadependency = Int[]
		i = 1
		for s in uniquespecies
			if haskey(isotope_mixtures, s)
				ind = findin(uniquespecies, isotope_mixtures[s])
				lind = length(ind)
				if lind == 1
					push!(deltaindex, i)
					push!(deltadependency, ind[1])
				elseif lind > 1
					error("More than one stable isotope dependency for $s")
					display(uniquespecies[ind])
					return
				end
			end
			i += 1
		end
		global deltaindex = deltaindex
		global deltadependency = deltadependency
	else
		warn("$filename is missing!")
	end
	filename = "data/cr-data-wells-$(probstamp).csv"
	if isfile(filename)
		rawwells = readcsv(filename)[2:end,[1,2,4,5,10]]
	else
		error("$filename is missing!")
		return
	end
	filename = "data/cr-data-pz-$(probstamp).csv"
	if isfile(filename)
		rawpz = readcsv(filename)[2:end,[1,2,4,5,10]]
	else
		error("$filename is missing!")
		return
	end
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
		if haskey(dict_species, longnames[i])
			names[i] = dict_species[longnames[i]]
		end
	end
	goodindices = 1:length(wells)
	if wellsset == ""
		# remove MCOI LAOI SCI R-6i TA-53i
		goodindices = filter(i->!contains(wells[i], "MCOI"), goodindices)
		goodindices = filter(i->!contains(wells[i], "LAOI"), goodindices)
		goodindices = filter(i->!contains(wells[i], "SCI"), goodindices)
		goodindices = filter(i->!contains(wells[i], "R-6i"), goodindices)
		goodindices = filter(i->!contains(wells[i], "TA-53i"), goodindices)
	else
		filename = "data/cr-wells-set$(wellsset).txt"
		if isfile(filename)
			ws = readdlm(filename)
			wu = unique(wells)
			sd = setdiff(ws, wu)
			if length(sd) > 0
				display(sd)
				error("There are wells in $filename missing!")
				return
			end
			sd = setdiff(wu, ws)
			if length(sd) > 0
				warn("The following wells will be removed!")
				display(sd)
			end
			for w in sd
				goodindices = filter(i->!contains(wells[i], w), goodindices)
			end
		else
			error("$filename is missing!")
			return
		end
	end
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
	name2j = Dict(zip(uniquespecies, 1:length(uniquespecies)))
	datacount = zeros(Int, length(uniquewells), length(uniquespecies))
	global datamatrix = zeros(Float32, length(uniquewells), length(uniquespecies))
	for index = 1:length(wells)
		i = wells2i[wells[index]]
		if haskey(name2j, names[index])
			j = name2j[names[index]]
			datacount[i, j] += 1
			datamatrix[i, j] += concs[index]
		end
	end
	wellorder, wellnameorder = getwellorder()
	info("Species ($(length(uniquespecies)))")
	display(uniquespecies)
	info("Wells ($(length(wellnameorder)))")
	display(wellnameorder)
	info("Observations count:")
	display([transposevector(["Wells"; uniquespecies]); wellnameorder datacount[wellorder,:]])
	info("Observations per well:")
	display([wellnameorder sum(datacount,2)[wellorder]])
	info("Observations per species:")
	display([uniquespecies sum(datacount,1)'])
	global datamatrix = convert(Array{Float32,2}, datamatrix ./ datacount) # gives NaN if there is no data, otherwise divides by the number of results
	info("Concentration matrix:")
	display([transposevector(["Wells"; uniquespecies]); wellnameorder datamatrix[wellorder,:]])
	info("Potential regularization penalty = $(sum(log(1.+abs(maximum(abs(datamatrix), 1))).^2))")
	global dataindex = collect(1:size(datamatrix, 2))
	global concindex = setdiff(dataindex, deltaindex)
	coord, coordheader = readdlm("data/cr-well-coord.dat", header=true)
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
	not_ok = false
	uniquelongnames = unique(longnames)
	for i in 1:length(uniquelongnames)
		if !haskey(dict_species, uniquelongnames[i])
			not_ok = true
			warn("Species name $(uniquelongnames[i]) in the data set is not defined in the dictionary!")
		end
	end
	if !not_ok
		println("ok")
	end
	
	info("Check species in the data set ...")
	not_ok = false
	for i in 1:length(uniquespecies_long)
		if indexin([uniquespecies_long[i]], uniquelongnames)[1] == 0
			not_ok = true
			warn("Species name $(uniquespecies_long[i]) defined in the dictionary is missing in the data set!")
		end
	end
	if !not_ok
		println("ok")
	end
end

"""
Perform rMF analyses
"""
function execute(range::Union{UnitRange{Int},Int}=1:maxbuckets; retries::Int=10, mixmatch::Bool=true, mixtures::Bool=true, normalize::Bool=false, scale::Bool=false, regularizationweight::Float32=convert(Float32, 0), weightinverse::Bool=false, matchwaterdeltas::Bool=false, quiet::Bool=true,clusterweights::Bool=true)
	if sizeof(datamatrix) == 0
		warn("Execute `rMF.loaddata()` first!")
		return
	end
	concmatrix = datamatrix
	if length(ratioindex) > 0
		concmatrix = datamatrix[:, concindex]
		ratiomatrix = datamatrix[:, ratioindex]
		nummixtures = size(concmatrix, 1)
		numconstituents = size(concmatrix, 2)
		info("Mixtures matrix:")
		display([transposevector(["Wells"; uniquespecies[concindex]]); uniquewells concmatrix])
		info("Ratio matrix:")
		display([transposevector(["Wells"; uniquespecies[ratioindex]]); uniquewells ratiomatrix])
		ratios = convert(Array{Float32, 3}, fill(NaN, nummixtures, numconstituents, numconstituents))
		for i = 1:nummixtures
			for j = 1:size(ratiocomponents, 2)
				a = ratiocomponents[1, j]
				b = ratiocomponents[2, j]
				ratios[i, a, b] = ratiomatrix[i]
			end
		end
	else
		ratios = nothing
	end
	if length(deltaindex) > 0
		concmatrix = datamatrix[:, concindex]
		deltamatrix = datamatrix[:, deltaindex]
		deltaindices = indexin(deltadependency, concindex)
		info("Mixtures matrix:")
		display([transposevector(["Wells"; uniquespecies[concindex]]); uniquewells concmatrix])
		info("Delta matrix:")
		display([transposevector(["Wells"; uniquespecies[deltaindex]]); uniquewells deltamatrix])
	else
		deltaindices = deltadependency
		deltamatrix = Array(Float32, 0, 0)
	end
	for numbuckets = range
		# NMFk using mixmatch 
		mixers[numbuckets], buckets[numbuckets], fitquality[numbuckets], robustness[numbuckets] = NMFk.execute(concmatrix, retries, numbuckets; deltas=deltamatrix, deltaindices=deltaindices, ratios=ratios, mixmatch=mixmatch, normalize=normalize, scale=scale, matchwaterdeltas=matchwaterdeltas, mixtures=mixtures, quiet=quiet, regularizationweight=regularizationweight, weightinverse=weightinverse, clusterweights=clusterweights)
		mixsum = sum(mixers[numbuckets], 2)
		index = find(mixsum .> 1.1) | find(mixsum .< 0.9)
		if length(index) > 0
			warn("The mixers do not add to 1")
			display(sum(mixers[numbuckets], 2))
			display(mixers[numbuckets])
		end
		println("Buckets = $numbuckets; Best objective function = $(fitquality[numbuckets]); Robustness = $(robustness[numbuckets])")
		if casekeyword == ""
			filename = "results/$case-$numbuckets-$retries.jld"
		else
			filename = "results/$case-$casekeyword-$numbuckets-$retries.jld"
		end
		JLD.save(filename, "wells", uniquewells, "species", uniquespecies, "mixers", mixers[numbuckets], "buckets", buckets[numbuckets], "fit", fitquality[numbuckets], "robustness", robustness[numbuckets], "regularizationweight", regularizationweight)
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
