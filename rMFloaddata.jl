"Check data for rMF analysis"
function check(casename::AbstractString, numruns::Int=100, keyword::AbstractString=""; noise::Bool=false, ns::Int=3, nw::Int=10, nc::Int=5, nd::Int=0, nr::Int=0, seed::Integer=0)
	if seed != 0
		srand(seed)
	end
	nb = max(ns, nc+nr+nd)
	numberofsourcesreconstruction = Array{Int64}(numruns)
	numberofsourcesrobustness = Array{Int64}(numruns)
	numberofsourcesaic = Array{Int64}(numruns)
	for i = 1:numruns
		info("$i")
		loaddata(casename, keyword; noise = noise, ns = ns, nw = nw, nc = nc, nd = nd, nr = nr, quiet = true)
		execute(1:nb, nooutput = true)
		# getresults(1:nb, brief = true)
		numberofsourcesreconstruction[i] = indmax(fitquality[1:nb-1] - fitquality[2:nb]) + 1
		numberofsourcesrobustness[i] = indmax(robustness[1:nb-1] - robustness[2:nb])
		numberofsourcesaic[i] = indmax(abs(aic[1:nb-1] - aic[2:nb])) + 1
	end
	cfit = count(i->numberofsourcesreconstruction[i].==ns, 1:numruns)
	crob = count(i->numberofsourcesrobustness[i].==ns, 1:numruns)
	caic = count(i->numberofsourcesaic[i].==ns, 1:numruns)
	info("Correct (Fit) = $cfit/$numruns")
	info("Correct (Silhouette) = $crob/$numruns")
	info("Correct (AIC) = $caic/$numruns")
	return
end

"Load data for rMF analysis"
function loaddata(casename::AbstractString, keyword::AbstractString=""; noise::Bool=false, ns::Int=3, nw::Int=10, nc::Int=5, nd::Int=0, nr::Int=0, quiet::Bool=false, seed::Integer=0)
	if seed != 0
		srand(seed)
	end
	global case = casename
	global casekeyword = keyword
	global mixers = Array{Array{Float64, 2}}(maxbuckets)
	global buckets = Array{Array{Float64, 2}}(maxbuckets)
	global fitquality = Array{Float64}(maxbuckets)
	global robustness = Array{Float64}(maxbuckets)
	global ratioindex = Int[]
	global deltaindex = Int[]
	global deltadependency = Int[]
	if casename == "test23delta"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["δA", "A", "B"]
		global datamatrix = convert(Array{Float32, 2}, [[0.1, 1.] [1., 0.1] [0.1, 1.]])
		global deltastandards = [1]
		global deltaindex = Int[1]
		global deltadependency = Int[2]
		global dataindex = collect(1:size(datamatrix, 2))
		global concindex = setdiff(dataindex, deltaindex)
		global wellcoord = [[0., 0.] [0., 100.]]
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test23delta2"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "δA"]
		global datamatrix = convert(Array{Float32, 2}, [[1., 0.1] [0.1, 1.] [0.1, 1.]])
		global deltaindex = Int[3]
		global deltadependency = Int[1]
		global dataindex = collect(1:size(datamatrix, 2))
		global concindex = setdiff(dataindex, deltaindex)
		global wellcoord = [[0., 0.] [0., 100.]]
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test23ratio"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "A/B"]
		global datamatrix = convert(Array{Float32, 2}, [[NaN, 2] [1, NaN] [1., 2.]])
		global ratioindex = Int[3]
		global ratiocomponents = Int[1, 2]
		global concindex = setdiff(collect(1:size(datamatrix,2)), ratioindex)
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.]]
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test27ratiosonly"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "C", "D", "A/B", "B/C", "C/D"]
		global datamatrix = convert(Array{Float32, 2}, [[NaN, NaN] [NaN, NaN] [NaN, NaN] [NaN, NaN] [1., 2.] [3., 2.] [4., 3.]])
		global ratioindex = Int[5,6,7]
		global ratiocomponents = Int[[1, 2] [2, 3] [3, 4]]
		global concindex = setdiff(collect(1:size(datamatrix,2)), ratioindex)
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.]]
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test23"
		global uniquewells = ["W1", "W2"]
		global uniquespecies = ["A", "B", "C"]
		global datamatrix = convert(Array{Float32,2}, [[1., 1.] [2., 4.] [2., 2.]])
		global concindex = collect(1:size(datamatrix,2))
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.]]
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test56s3"
		global uniquewells = ["W1", "W2", "W3", "W4", "W5"]
		global uniquespecies = ["A", "B", "C", "D", "E", "F"]
		global truemixer = [ 0.1 0.3 0.7;
							 0.7 0.3 0.1;
							 0.3 0.3 0.4;
							 0.0 0.0 1.0;
							 1.0 0.0 0.0]
		global truebucket = [ 1.0 0.0 0.0 0.1 0.1 0.1;
							  0.0 1.0 0.0 0.1 0.1 0.1;
							  0.0 0.0 1.0 0.1 0.1 0.1]
		if noise
			noise_matrix = randn(length(uniquewells), length(uniquespecies)) / 100
		else
			noise_matrix = 0
		end
		global datamatrix = convert(Array{Float32,2}, truemixer * truebucket + noise_matrix)
		global concindex = collect(1:size(datamatrix,2))
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.] [50., 50.] [0., 50.] [50., 0.]]
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test56s4"
		global uniquewells = ["W1", "W2", "W3", "W4", "W5"]
		global uniquespecies = ["A", "B", "C", "D", "E", "F"]
		global truemixer = [ 0.1 0.7 0.1 0.1;
							 0.2 0.4 0.1 0.3;
							 0.1 0.1 0.7 0.1;
							 0.4 0.2 0.1 0.3;
							 0.7 0.1 0.1 0.1]
		global truebucket = [ 100.0   0.1   0.1   0.1  10.0   0.0;
								0.1 100.0   0.3   0.4  90.0  30.0;
								0.5   0.4 100.0   0.2  90.0   0.0;
								0.1   0.4   0.3 100.0  10.0  30.0]
		if noise
			noise_matrix = randn(length(uniquewells), length(uniquespecies)) / 100
		else
			noise_matrix = 0
		end
		global datamatrix = convert(Array{Float32,2}, truemixer * truebucket + noise_matrix)
		global concindex = collect(1:size(datamatrix,2))
		global dataindex = concindex
		global wellcoord = [[0., 0.] [0., 100.] [50., 50.] [0., 50.] [50., 0.]]
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	if casename == "test"
		if ns < 1
			error("Number of sources should be larger than zero!")
		end
		if nw > 0
			global uniquewells = map(i->"W$i", 1:nw)
		else
			error("Number of wells should be larger than zero!")
			return
		end
		if nc > 0
			global uniquespecies = map(i->string(Char(65 + (i-1)%26))^Int(ceil(i/26)), 1:nc)
		else
			error("Number of wells should be larger than zero!")
			return
		end
		global concindex = collect(1:nc)
		global dataindex = collect(1:(nc+nd))
		global wellcoord = []
		mixer = rand(nw, ns)
		mixer = (mixer .* 2).^10
		mixer_norm = diagm(1 ./ vec(sum(mixer, 2)))
		global truemixer = mixer_norm * mixer
		global truemixer = truemixer ./ (sum(truemixer, 2))
		bucket = rand(ns, nc)
		bucket = (bucket .* 2).^4
		bucket_norm = diagm(1 ./ vec(maximum(bucket, 1)))
		global truebucket = bucket * bucket_norm
		if noise
			noise_matrix = randn(nw, nc + nd + nr) / 100
		else
			noise_matrix = 0
		end
		if nd > 0
			if nd > nc
				warn("Number of deltas cannot be larger than the number of concentrations!")
				return
			end
			deltas = map(i->"δ" * string(Char(65 + (i-1)%26))^Int(ceil(i/26)), 1:nd)
			global uniquespecies = vcat(uniquespecies, deltas)
			global deltaindex = map(i->(nc + i), 1:nd)
			global deltadependency = map(i->i, 1:nd)
			global deltastandards = map(i->1, 1:nd)
			# deltas = (rand(ns, nd) .* 2).^10
			deltas = rand(ns, nd)
			deltas_norm = diagm(1 ./ vec(maximum(deltas, 1)))
			global truedeltas = deltas * deltas_norm
			deltas = NMFk.computedeltas(truemixer, truebucket, truedeltas, indexin(deltadependency, concindex))
			global datamatrix = convert(Array{Float32,2}, hcat(truemixer * truebucket, deltas) + noise_matrix)
			datamatrix[:,1:nd] = NaN
		end
		if nr > 0
			if nr > nc - 1
				warn("Number of deltas cannot be larger than the number of concentrations minus one!")
				return
			end
			ratios = map(i->string(Char(65 + (i-1)%26))^Int(ceil(i/26)) * "/" * string(Char(65 + (i)%26))^Int(ceil(i/26)), 1:nr)
			global uniquespecies = vcat(uniquespecies, ratios)
			global ratioindex = map(i->(nc + i), 1:nr)
			global ratiocomponents = [1:nr 2:nr+1]'
			global datamatrix = convert(Array{Float32,2}, truemixer * truebucket)
			global trueratios = map(Float32, (datamatrix[i,j] / datamatrix[i,j + 1]) for i=1:nw, j=1:nr)
			global datamatrix = convert(Array{Float32,2}, hcat(datamatrix, trueratios) + noise_matrix)
			datamatrix[:,1:nr+1] = NaN
		end
		if nd <= 0 && nr <=0
			global datamatrix = convert(Array{Float32,2}, truemixer * truebucket + noise_matrix)
		end
		global case = casename * "_" * string(nw) * "_" * string(nc) * "_" * string(nd) * "_" * string(nr) * "_" * string(ns)
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	end
	filename = "data/" * casename * ".csv"
	if isfile(filename)
		rawdata = readcsv(filename)
		rawdata[rawdata .== " "] = NaN
		rawdata[rawdata .== ""] = NaN
		global uniquewells = rawdata[2:end, 1]
		global uniquespecies = rawdata[1, 2:end]
		global datamatrix = convert(Array{Float32,2}, rawdata[2:end, 2:end])
		global concindex = collect(1:size(datamatrix,2))
		global dataindex = concindex
		global truebucket = Array{Float64}(0)
		global truedeltas = Array{Float64}(0)
		global trueratios = Array{Float64}(0)
		global truemixer = Array{Float64}(0)
		info("Species ($(length(uniquespecies)))")
		display(uniquespecies)
		info("Wells ($(length(uniquewells)))")
		display(uniquewells)
		!quiet && info("Concentration matrix:")
		!quiet && display([transposevector(["Wells"; uniquespecies]); uniquewells datamatrix])
		return
	else
		warn("File $filename is missing!")
	end
end

"Get well order"
function getwellorder(sort::Bool=false)
	nw = length(uniquewells)
	if nw == 0
		return
	end
	if !sort
		if isfile("data/cr-well-order-WE.dat")
			wells2i = Dict(zip(uniquewells, 1:nw))
			wellnameorder = strip.(readdlm("data/cr-well-order-WE.dat"))
			nwellnameorder = length(wellnameorder)
			info("Number of wells in data/cr-well-order-WE.dat is $nwellnameorder")
			wellorder = zeros(Int, length(wellnameorder))
			countmissing = 0
			for i = 1:length(wellorder)
				if haskey(wells2i, wellnameorder[i])
					wellorder[i] = wells2i[wellnameorder[i]]
				else
					println("Well $(wellnameorder[i]) is missing in the input data set")
					countmissing += 1
				end
			end
			info("Number of wells in data/cr-well-order-WE.dat missing in the input data set is $countmissing")
			wellmissing = wellorder .== 0
			indexmissing = find(wellmissing)
			wellavailable = wellorder .!= 0
			indexavailale = find(wellavailable)
			wellorder = wellorder[indexavailale]
			wellnameorder = wellnameorder[indexavailale]
		else
			warn("data/cr-well-order-WE.dat is missing!")
			wellorder = collect(1:nw)
			wellnameorder = uniquewells
		end
		if nw > 0 && length(wellnameorder) > 0 && nw != length(wellnameorder)
			if nw > nwellnameorder
				wellorder = sortperm(uniquewells)
				wellnameorder = uniquewells[wellorder]
				warn("There are wells missing in data/cr-well-order-WE.dat")
				warn("Original order preserved!")
			end
		end
	else
		wellorder = sortperm(uniquewells)
		wellnameorder = uniquewells[wellorder]
	end
	return wellorder, wellnameorder
end

"Display rMF data"
function displayconc(name::AbstractString)
	displayconc([name])
end
function displayconc(names::Vector{String})
	if length(names) == 0
		display([transposevector(["Wells"; uniquespecies]); wellnameorder datamatrix[wellorder,:]])
	else
		i = findin(uniquespecies, names)
		if length(i) > 0
			display([transposevector(["Wells"; uniquespecies[i]]); wellnameorder datamatrix[wellorder,i]])
		else
			i = findin(wellnameorder, names)
			if length(i) > 0
				display(transposematrix([transposevector(["Wells"; uniquespecies]); wellnameorder[i] datamatrix[wellorder[i],:]]))
			end
		end
	end
end

"Load data for rMF analysis"
function loaddata(probstamp::Int64=20160102, keyword::AbstractString=""; wellsset::AbstractString="", speciesset::AbstractString="", datemin::Date=Date(1964), datemax::Date=Date(2064), quiet::Bool=false)
	casestring = keyword
	global truebucket = Array{Float64}(0)
	global truedeltas = Array{Float64}(0)
	global trueratios = Array{Float64}(0)
	global truemixer = Array{Float64}(0)
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
	global mixers = Array{Array{Float64, 2}}(maxbuckets)
	global buckets = Array{Array{Float64, 2}}(maxbuckets)
	global fitquality = Array{Float64}(maxbuckets)
	global robustness = Array{Float64}(maxbuckets)
	global ratioindex = Int[]

	# read well data
	filename = "data/cr-data-wells-$(probstamp).csv"
	if isfile(filename)
		rawwells = readcsv(filename)[2:end,[1,2,4,5,10]]
	else
		error("$filename is missing!")
		return
	end

	# read pz data
	filename = "data/cr-data-pz-$(probstamp).csv"
	if isfile(filename)
		rawpz = readcsv(filename)[2:end,[1,2,4,5,10]]
	else
		error("$filename is missing!")
		return
	end

	# process the data
	rawdata = [rawwells; rawpz]
	wells = rawdata[:, 1]
	dates = rawdata[:, 2]
	longnames = rawdata[:, 3]
	names = rawdata[:, 4]
	concs = rawdata[:, 5]
	@assert length(longnames) == length(names)
	@assert length(names) == length(concs)
	info("Total data record count $(length(longnames))")
	inputspecies = sort(unique(longnames))
	info("Species in the input file $filename:")
	display(inputspecies)

	# read species renaming dictionary
	filename = "data/cr-species.jld"
	if isfile(filename)
		global dict_species = JLD.load(filename, "species")
	else
		error("$filename is missing!")
		return
	end

	# set the "cool" (short) species names
	for i = 1:length(longnames)
		if haskey(dict_species, longnames[i])
			names[i] = dict_species[longnames[i]]
		end
	end
	global uniquespecies = unique(names)

	# sort/filter species based on an input species set file
	if speciesset != ""
		filename = "data/cr-species-set$(speciesset).txt"
		if isfile(filename)
			ss = readdlm(filename)
			sd = setdiff(ss, uniquespecies)
			if length(sd) > 0
				warn("There are species in $filename that are not defined in the input data set!")
				display(sd)
			end
			sd = setdiff(uniquespecies, ss)
			if length(sd) > 0
				warn("The following species will be ignored if provided in the data set!")
				display(sd)
			end
			ind = findin(ss, uniquespecies)
			global uniquespecies = ss[ind]
			dnames = ss
		else
			error("$filename is missing!")
		end
	end

	filename = "data/cr-stable-isotope-mixtures.jld"
	if isfile(filename)
		deltamixtures = JLD.load(filename, "deltamixtures")
		deltaindex = Int[]
		deltadependency = Int[]
		i = 1
		for s in uniquespecies
			if haskey(deltamixtures, s)
				ind = findin(uniquespecies, deltamixtures[s])
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
	filename = "data/cr-stable-isotope-standards.jld"
	if isfile(filename)
		dict_deltastandards = JLD.load(filename, "deltastandards")
		deltastandards = Float64[]
		for s in uniquespecies[deltaindex]
			if haskey(dict_deltastandards, s)
				push!(deltastandards, dict_deltastandards[s])
			else
				throw("Delta Standard is missing for $(s)")
			end
		end
		global deltastandards = deltastandards
	else
		warn("$filename is missing!")
	end

	# rename chromium wells
	goodindices = 1:length(wells)
	for i in goodindices
		m = match(r"([A-Z].*) S([1-9])", wells[i])
		if m != nothing && length(m.captures) == 2
			wells[i] = m.captures[1] * "_" * m.captures[2]
		end
		m = match(r"(C[rR][Ee][Xx])(.*)", wells[i])
		if m != nothing && length(m.captures) == 2
			wells[i] = "Ex" * m.captures[2]
		end
		m = match(r"(C[rR][Ii][Nn])(.*)", wells[i])
		if m != nothing && length(m.captures) == 2
			wells[i] = "In" * m.captures[2]
		end
		m = match(r"(C[rR][Pp][Zz])(.*)", wells[i])
		if m != nothing && length(m.captures) == 2
			wells[i] = "Pz" * m.captures[2]
		end
	end
	if wellsset == ""
		# remove MCOI LAOI SCI R-6i TA-53i
		sd = ["MCOI", "LAOI", "SCI", "R-6i", "TA-53i"]
		warn("The following wells will be removed!")
		display(sd)
		for w in sd
			goodindices = filter(i->(wells[i] != w), goodindices)
		end
	else
		filename = "data/cr-wells-set$(wellsset).txt"
		if isfile(filename)
			ws = readdlm(filename)
			wu = unique(wells)
			sd = setdiff(ws, wu)
			if length(sd) > 0
				info("Wells in the data set:")
				display(wu)
				info("Wells in $filename:")
				display(ws)
				if length(sd) > 0
					warn("There are wells in the input data set that are missing in $(filename)!")
					info("Missing wells in $(filename):")
					display(sd)
				end
			end
			sd = setdiff(wu, ws)
			if length(sd) > 0
				warn("The following wells will be removed!")
				display(sd)
			end
			for w in sd
				goodindices = filter(i->(wells[i] != w), goodindices)
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
	datetime=Array{Date}(length(dates))
	for i = 1:length(dates)
		d = 0
		try
			d = Dates.Date(dates[i], "mm/dd/yyyy")
		catch
			try
				d = Dates.Date(dates[i], "mm-dd-yyyy")
			catch
				@show dates[i]
			end
		end
		datetime[i] = d
	end
	goodindices = filter(i->(datetime[i] >= datemin && datetime[i] <= datemax), goodindices)
	dates = dates[goodindices]
	wells = wells[goodindices]
	names = names[goodindices]
	longnames = longnames[goodindices]
	concs = concs[goodindices]
	info("Number of processed data entries: $(length(concs))")

	global uniquewells = unique(wells)
	wells2i = Dict(zip(uniquewells, 1:length(uniquewells)))
	name2j = Dict(zip(uniquespecies, 1:length(uniquespecies)))
	datacount = zeros(Int, length(uniquewells), length(uniquespecies))
	global datamatrix = zeros(Float32, length(uniquewells), length(uniquespecies))
	global sdmatrix = zeros(Float32, length(uniquewells), length(uniquespecies))
	for index = 1:length(wells)
		i = wells2i[wells[index]]
		if haskey(name2j, names[index])
			j = name2j[names[index]]
			datacount[i, j] += 1
			datamatrix[i, j] += concs[index]
			sdmatrix[i, j] += concs[index] * concs[index]
		end
	end

	wellorder, wellnameorder = getwellorder()
	global wellorder
	global wellnameorder

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

	sdmatrix = sqrt(abs(sdmatrix - (datamatrix .^2) ./ datacount))
	sdmatrix[isnan(sdmatrix)] = 0
	sdmatrix[datacount.==1] = 0
	global datamatrix = convert(Array{Float32,2}, datamatrix ./ datacount) # gives NaN if there is no data, otherwise divides by the number of results

	!quiet && info("Concentration matrix:")
	!quiet && display([transposevector(["Wells"; uniquespecies]); wellnameorder datamatrix[wellorder,:]])
	!quiet && info("Concentration standard deviation matrix:")
	!quiet && display([transposevector(["Wells"; uniquespecies]); wellnameorder sdmatrix[wellorder,:]])

	!quiet && info("Maximum standard deviation for various species:")
	!quiet && display([uniquespecies maximum(sdmatrix, 1)'])
	sdmatrix2 = deepcopy(sdmatrix)
	info("Largest standard deviations:")
	for i = 1:20
		indmaxsd = ind2sub(size(sdmatrix2), indmax(sdmatrix2))
		mmm = sdmatrix2[indmaxsd[1],indmaxsd[2]]
		println("$i - $(uniquewells[indmaxsd[1]]) / $(uniquespecies[indmaxsd[2]]): standard deviations $(mmm) sample size $(datacount[indmaxsd[1],indmaxsd[2]])")
		if mmm <= 0
			break
		end
		sdmatrix2[indmaxsd[1],indmaxsd[2]] = 0
	end
	sdmatrix2 = deepcopy(sdmatrix)
	info("Smallest non-zero standard deviations:")
	for i = 1:5
		indminsd = ind2sub(size(sdmatrix2), indmin(sdmatrix2[sdmatrix2.>0]))
		mmm = sdmatrix2[indminsd[1],indminsd[2]]
		println("$i - $(uniquewells[indminsd[1]]) / $(uniquespecies[indminsd[2]]): standard deviations $(mmm) sample size $(datacount[indminsd[1],indminsd[2]])")
		sdmatrix2[indminsd[1],indminsd[2]] = Inf
	end
	info("Potential regularization penalty = $(sum(log(1 .+ abs(maximum(datamatrix, 1))).^2))")
	global dataindex = collect(1:size(datamatrix, 2))
	global concindex = setdiff(dataindex, deltaindex)
	coord, coordheader = readdlm("data/cr-well-coord.dat", header=true)
	global wellcoord = Array{Float64}(length(uniquewells), 2)
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
	for i = 1:length(uniquelongnames)
		if !haskey(dict_species, uniquelongnames[i])
			not_ok = true
			warn("Species name `$(uniquelongnames[i])` in the data set is not defined in the dictionary!")
		end
	end
	if !not_ok
		println("ok")
	end

	info("Check species in the data set ...")
	not_ok = false
	uniquespecies_long = collect(keys(dict_species))
	for i = 1:length(uniquespecies_long)
		if indexin([uniquespecies_long[i]], uniquelongnames)[1] == 0
			not_ok = true
			warn("Species name `$(uniquespecies_long[i])` defined in the dictionary is missing in the input data set!")
		end
	end
	if !not_ok
		println("ok")
	end

	info("Check species for undefined species the data set ...")
	not_ok = false
	speciescount = sum(datacount,1)'
	badspeciesindex = speciescount .== 0
	badspecies = uniquespecies[badspeciesindex]
	if length(badspecies) > 0
		warn("Species should be removed because they have no data")
		display(badspecies)
		not_ok = true
	end
	if !not_ok
		println("ok")
	end
end
