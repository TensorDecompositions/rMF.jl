"Retrieve saved rMF results"
function getresults(range::Union{UnitRange{Int},Int}=1:maxbuckets, keyword::AbstractString=""; retries::Int=10, brief::Bool=false)
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

	numwells = size(datamatrix, 1)
	numconstituents = size(datamatrix, 2)
	@assert numwells == length(wellnameorder)
	@assert numconstituents == length(uniquespecies)
	for numbuckets in range
		filename = "results/$(casestring)-$numbuckets-$retries.jld"
		if isfile(filename)
			j = JLD.load(filename)
			mixers[numbuckets] = j["mixers"]
			buckets[numbuckets] = j["buckets"]
			fitquality[numbuckets] = j["fit"]
			robustness[numbuckets] = j["robustness"]
			if haskey(j, "aic")
				aic[numbuckets] = j["aic"]
			else
				aic[numbuckets] = NaN
			end
			order = DataStructures.OrderedDict(zip(uniquespecies, 1:length(uniquespecies)))
			if haskey(j, "uniquewells") && haskey(j, "uniquespecies")
				wells = j["uniquewells"]
				species = j["uniquespecies"]
				remap = DataStructures.OrderedDict(zip(species, 1:length(species)))
			else
				if isfile("data/cr-species-order-$(case).jld")
					dictold=JLD.load("data/cr-species-order-$(case).jld","dictionary")
					remap = DataStructures.OrderedDict(zip(collect(keys(dictold)), 1:length(uniquespecies)))
				else
					remap = DataStructures.OrderedDict(zip(uniquespecies, 1:length(uniquespecies)))
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

		orderedbuckets = similar(buckets[numbuckets])
		global spredictions = Array{Array{Float64, 2}}(numbuckets)
		if length(deltaindex) > 0
			numdeltas = length(deltaindex)
			numconc = numconstituents - numdeltas
			H_conc = buckets[numbuckets][:,1:numconc]
			H_deltas = buckets[numbuckets][:,numconc+1:end]
			orderedbuckets[:, concindex] = H_conc
			orderedbuckets[:, deltaindex] = H_deltas
			global predictions = similar(datamatrix)
			predictions[:, concindex] = mixers[numbuckets] * H_conc
			predictions[:, deltaindex] = NMFk.computedeltas(mixers[numbuckets], H_conc, H_deltas, indexin(deltadependency, concindex))
			tpredictions = zeros(size(datamatrix))
			for i = 1:numbuckets
				spredictions[i] = similar(datamatrix)
				spredictions[i][:, concindex] = mixers[numbuckets][:,i:i] * H_conc[i:i,:]
				spredictions[i][:, deltaindex] = NMFk.computedeltas(reshape(mixers[numbuckets][:,i:i], numwells, 1), H_conc[i:i,:], H_deltas[i:i,:], indexin(deltadependency, concindex), compute_contributions=true)
				tpredictions = tpredictions + spredictions[i]
			end
			tpredictions[:, deltaindex] .*= predictions[:, deltaindex] ./ tpredictions[:, deltaindex]
			@Base.Test.test isapprox(maximum(abs(predictions .- tpredictions)), 0, atol=1e-3)
		else
			predictions = mixers[numbuckets] * buckets[numbuckets]
			orderedbuckets = buckets[numbuckets]
			for i = 1:numbuckets
				spredictions[i] = similar(datamatrix)
				spredictions[i] = mixers[numbuckets][:,i:i] * buckets[numbuckets][i:i,:]
			end
		end
		if length(ratioindex) > 0 # Compute ratios
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
			for k = 1:numbuckets
				p = similar(datamatrix)
				p[:, concindex] = spredictions[k]
				for i = 1:numwells
					for j = 1:size(ratiocomponents, 2)
						a = ratiocomponents[1, j]
						b = ratiocomponents[2, j]
						p[i, ratioindex[j]] = p[i, a] / p[i, b]
					end
				end
				spredictions[k] = p
			end
		end

		indexnan = isnan.(datamatrix)
		d = copy(datamatrix)
		d[indexnan] = 0
		try
			@assert size(d) == size(predictions)
		catch
			error("There is a model setup mismatch! Data and prediction matrices do not match!\nData matrix: $(size(d))\nPrediction matrix: $(size(predictions))")
		end
		errors = d - predictions
		errors[indexnan] = 0.
		of = sum(errors.^2)
		errors[indexnan] = NaN
		relerrors = errors ./ d
		relerrors[indexnan] = NaN
		regularization_penalty = sum(log.(1.+abs.(orderedbuckets)).^2) / numbuckets

		vector_errors = vec(errors[.!indexnan])
		stddeverrors = std(vector_errors)
		if stddeverrors > 0
			kstest = HypothesisTests.ExactOneSampleKSTest(vector_errors, Distributions.Normal(0, stddeverrors))
			kstestdelta = kstest.δ
			kstestpval = HypothesisTests.pvalue(kstest)
			if kstestdelta > kstestpval
				kstestvertdict = "not normal"
			else
				kstestvertdict = "normal"
			end
		else
			kstestdelta = 0
			kstestpval = Inf
			kstestvertdict = "unknown"
		end

		f = open("results/$(casestring)-$numbuckets-$retries-stats.dat", "w")
		println(f, "Number of buckets: $(numbuckets)")
		println(f, "* Fit: $(fitquality[numbuckets])")
		println(f, "* Fit check: $(of)")
		println(f, "* Silhouette: $(robustness[numbuckets])")
		println(f, "* KS Test: $kstestvertdict ($(kstestdelta) $(kstestpval))")
		println(f, "* AIC: $(aic[numbuckets])")
		println(f, "* Error stats:")
		println(f, "  - Mean: $(mean(vector_errors))")
		println(f, "  - Variance: $(var(vector_errors))")
		println(f, "  - Standard Deviation: $(stddeverrors)")
		println(f, "  - Maximum: $(maximum(vector_errors))")
		println(f, "  - Minimum: $(minimum(vector_errors))")
		println(f, "  - Skewness: $(StatsBase.skewness(vector_errors))")
		println(f, "  - Kurtosis: $(StatsBase.kurtosis(vector_errors))")
		close(f)

		if brief
			print("Sources: $(@sprintf("%2d", numbuckets)) Fit: $(@sprintf("%12.7g", fitquality[numbuckets])) ")
			if of - fitquality[numbuckets] > 1e-1
				print("(check fails: $(@sprintf("%12.7g", of))) ")
			end
			println("Silhouette: $(@sprintf("%12.7g", robustness[numbuckets])) AIC: $(@sprintf("%12.7g", aic[numbuckets])) KS: $(@sprintf("%12.7g", kstestdelta)) StdDev: $(@sprintf("%12.7g", stddeverrors))")
			continue
		else
			info("Fit quality: $(fitquality[numbuckets]) (check = $(of)) (regularization penalty = $(regularization_penalty))")
			if of - fitquality[numbuckets] > 1e-1
				warn("Objective function test fails!")
			end
			info("Silhouette: $(robustness[numbuckets])")
		end

		info("Match error statistics:")
		println("KS Test: $kstestvertdict ($(kstestdelta) $(kstestpval))")
		println("AIC: $(aic[numbuckets])")
		println("Mean: $(mean(vector_errors))")
		println("Variance: $(var(vector_errors))")
		println("Standard Deviation: $(stddeverrors)")
		println("Maximum: $(maximum(vector_errors))")
		println("Minimum: $(minimum(vector_errors))")
		println("Skewness: $(StatsBase.skewness(vector_errors))")
		println("Kurtosis: $(StatsBase.kurtosis(vector_errors))")

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
		MA = Array{Compose.Context}(MArows, MAcols)
		i = 1
		for w in wellorder
			b = abs.(hcat(map(i->collect(spredictions[i][w,:]), 1:numbuckets)...)) ./ abs.(predictions[w:w, :]')
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
		Gadfly.draw(Gadfly.SVG(filename, MArows * (0.8Gadfly.inch + numbuckets * 0.2Gadfly.inch), MAcols * (0.4Gadfly.inch + numconstituents * 0.2Gadfly.inch)), gs)
		filename = "results/$(casestring)-$numbuckets-$retries-wellmixtures.png"
		Gadfly.draw(Gadfly.PNG(filename, MArows * (0.8Gadfly.inch + numbuckets * 0.2Gadfly.inch), MAcols * (0.4Gadfly.inch + numconstituents * 0.2Gadfly.inch)), gs)

		info("Match errors:")
		display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder errors[wellorder, :]]))
		f = open("results/$(casestring)-$numbuckets-$retries-errors.dat", "w")
		writedlm(f, transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder errors[wellorder, :]]))
		close(f)

		info("Histogram of the estimation errors:")
		g = Gadfly.plot(x=vector_errors, Gadfly.Geom.histogram())
		filename = "results/$(casestring)-$numbuckets-$retries-error_histogram.png"
		Gadfly.draw(Gadfly.PNG(filename, 6Gadfly.inch, 4Gadfly.inch), g)
		if isdefined(:madsdisplay)
			madsdisplay(filename)
		end

		indmaxerror = ind2sub(size(errors), indmax(abs.(errors)))
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

		indmaxerror = ind2sub(size(relerrors), indmax(abs.(relerrors)))
		info("The largest absolute relative match error is for $(wellnameorder[indmaxerror[1]]) / $(uniquespecies[indmaxerror[2]]).")
		println("Observation: $(datamatrix[indmaxerror...])")
		println("Prediction: $(predictions[indmaxerror...])")
		println("Error: $(errors[indmaxerror...])")
		println("Relative error: $(relerrors[indmaxerror...])")
		display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder[indmaxerror[1]] relerrors[indmaxerror[1]:indmaxerror[1], :]]))

		info("Max/min Species in Sources per species:")
		maxs1 = maximum(orderedbuckets, 1)
		mins1 = minimum(orderedbuckets, 1)
		display([uniquespecies[dataindex] maxs1' mins1'])

		info("Max/min Species in Sources per buckets:")
		maxs2 = maximum(orderedbuckets, 2)
		mins2 = minimum(orderedbuckets, 2)
		display([maxs2 mins2])

		info("Mixers:")
		display([uniquewells mixers[numbuckets]])

		info("Sources:")
		display([uniquespecies[dataindex] orderedbuckets[source_index,:]'])
		f = open("results/$(casestring)-$numbuckets-$retries-buckets.dat", "w")
		writedlm(f, [uniquespecies[dataindex] orderedbuckets[source_index,:]'])
		close(f)

		info("Normalized buckets:")
		for i = 1:length(maxs1)
			if maxs1[i] == mins1[i]
				mins1[i] = 0
			end
		end
		for i = 1:length(maxs2)
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
		bucketimpact = Array{Float64}(numbuckets, numconstituents)
		for s = 1:numconstituents
			for i = 1:numbuckets
				# bucketimpact[i, s] = sum(abs((spredictions[i][:, s])./predictions[:,s]))
				bucketimpact[i, s] = sum(abs.((spredictions[i][:, s])))
			end
		end
		bucketimpactwells = Array{Float64}(numbuckets, numwells)
		for w = 1:numwells
			for i = 1:numbuckets
				# bucketimpactwells[i, w] = sum(abs((spredictions[i][w, :])./predictions[w,:]))
				bucketimpactwells[i, w] = sum(abs.((spredictions[i][w, :])))
			end
		end

		info("Ordered buckets to capture the overall impact on the species concentrations:")
		display([uniquespecies[dataindex] bucketimpact[source_index, dataindex]'])

		info("Max/min Species in model predictions for each bucket:")
		maxm = maximum(bucketimpact, 1)
		minm = minimum(bucketimpact, 1)
		display([uniquespecies maxm' minm'])
		for i = 1:length(maxm)
			if maxm[i] == minm[i]
				minm[i] = 0
			end
		end

		info("Max/min Species in model predictions for each well:")
		maxiw = maximum(bucketimpactwells, 1)
		miniw = minimum(bucketimpactwells, 1)
		display([wellnameorder maxiw' miniw'])
		for i = 1:length(maxiw)
			if maxiw[i] == miniw[i]
				miniw[i] = 0
			end
		end

		bucketimpactwells[source_index, wellorder] = (bucketimpactwells .- miniw) ./ (maxiw - miniw)
		gmixers = Gadfly.spy(bucketimpactwells', Gadfly.Scale.y_discrete(labels = i->wellnameorder[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Wells"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpactwells.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpactwells.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)

		info("Ordered buckets normalized to capture the overall impact on the species concentrations:")
		bucketimpact[source_index, :] = (bucketimpact .- minm) ./ (maxm - minm)
		display([uniquespecies bucketimpact'])

		gbucket = Gadfly.spy(bucketimpact', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpact.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
		filename = "results/$(casestring)-$numbuckets-$retries-bucketimpact.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)

		# s1buckets[s1buckets.<1e-6] = 1e-6
		gbucket = Gadfly.spy(s1buckets', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)

	#=
		gbucket = Gadfly.spy(s2buckets', Gadfly.Scale.y_discrete(labels = i->uniquespecies[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Species"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets2.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
		filename = "results/$(casestring)-$numbuckets-$retries-buckets2.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gbucket)
	=#
		info("Ordered buckets:")
		display([uniquespecies[dataindex] orderedbuckets[source_index,:]'])

		if sizeof(truebucket) > 0
			info("True buckets:")
			display([uniquespecies[dataindex] truebucket'])
		end

		info("Ordered mixers:")
		smixers = mixers[numbuckets][wellorder, source_index]
		smixers[smixers .< 0] = 0
		smixers[smixers .> 1] = 1
		display([wellnameorder smixers])

		if sizeof(truemixer) > 0
			info("True mixers:")
			display([wellnameorder truemixer])
		end

		gmixers = Gadfly.spy(smixers, Gadfly.Scale.y_discrete(labels = i->wellnameorder[i]), Gadfly.Scale.x_discrete,
					Gadfly.Guide.YLabel("Wells"), Gadfly.Guide.XLabel("Sources"), Gadfly.Guide.colorkey(""),
					Gadfly.Theme(point_size=20Gadfly.pt, major_label_font_size=14Gadfly.pt, minor_label_font_size=12Gadfly.pt, key_title_font_size=16Gadfly.pt, key_label_font_size=12Gadfly.pt),
					Gadfly.Scale.ContinuousColorScale(Gadfly.Scale.lab_gradient(parse(Colors.Colorant, "green"), parse(Colors.Colorant, "yellow"), parse(Colors.Colorant, "red")), minvalue = 0, maxvalue = 1))
		# filename, format = Mads.setimagefileformat(filename, format)
		filename = "results/$(casestring)-$numbuckets-$retries-mixers.svg"
		Gadfly.draw(Gadfly.SVG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)
		filename = "results/$(casestring)-$numbuckets-$retries-mixers.png"
		Gadfly.draw(Gadfly.PNG(filename,6Gadfly.inch,12Gadfly.inch), gmixers)

		if !isdefined(rMF, :wellcoord) || length(wellcoord) == 0
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
		zi = Array{Float64}(length(xi), length(yi))
		wm = mixers[numbuckets][goodindices,source_index]
		for s = 1:numbuckets
			z = wm[:, s]
			z[z.<1e-6] = 1e-6
			z = log10.(z)
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
			PyPlot.scatter(wc[:,1], wc[:,2], marker="o", c=z, s=70, cmap="jet")
			PyPlot.scatter(wc[:,1], wc[:,2], marker="o", c=z, s=68, cmap="jet")
			PyPlot.clim(zmin, zmax)
			for i = 1:length(wn)
				PyPlot.annotate(wn[i], xy=(wc[i,1], wc[i,2]), xytext=(-2, 2), fontsize=8, color="white", textcoords="offset points", ha="right", va="bottom", path_effects=[PathEffects.withStroke(linewidth=1, foreground="black")])
			end
			PyPlot.yticks([538500, 539500], ["538500", "539500"], rotation="vertical")
			PyPlot.tight_layout()
			PyPlot.savefig("results/$(casestring)-$numbuckets-$retries-source-$s.png")
			PyPlot.close()
		end
		#=
		selectedspecies =[1,3,5,7]
		swm = map(i->rMF.mixers[5][goodindices,i] * rMF.buckets[5][i,selectedspecies], source_index)
		for b = 1:numbuckets
			for s = 1:length(selectedspecies)
				z in swm[b][:, s]
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
		=#
	end
end