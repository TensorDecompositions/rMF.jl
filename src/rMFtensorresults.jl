function plottensorresults(X::Array, Xe::Array, W::Array, ns::Number=size(W)[2]; figuredir="tensor-figures", isn=isnan.(X), keyword="", prefix="$(rMF.case)-w$(rMF.wellssetid)-s$(rMF.speciessetid)-y$(rMF.tensoryearstep)-$keyword$ns-")
	sourcenames = map(i->"S$i", 1:ns)

	Z = deepcopy(X)
	Ze = deepcopy(Xe)
	Z[isn] .= 0
	Ze[isn] .= 0
	# [vec(maximum(maximum(Z, 3), 1)) vec(maximum(maximum(Ze, 3), 1))]
	smax = vec(maximum(Z, (1,3)))
	smaxw = max.(maximum(Z, 3), maximum(Ze, 3))
	Z[isn] .= Inf
	Ze[isn] .= Inf
	# [vec(minimum(minimum(Z, 3), 1)) vec(minimum(minimum(Ze, 3), 1))]
	smin = vec(minimum(Z, (1,3)))
	sminw = min.(minimum(Z, 3), minimum(Ze, 3))
	E = abs.(X[.!isn] .- Xe[.!isn])
	info("Maximum error: $(maximum(E))")
	Xn = deepcopy(X)
	for i=1:length(smax)
		Xn[:, i, :] ./= smax[i]
	end

	scolors = ["red", "blue", "green", "orange", "purple", "brown", "cyan", "magenta"]

	info("Plot species concentrations for each well ...")
	for (c, w) in enumerate(rMF.wellnameorder)
		p = Vector{Any}(length(rMF.uniquespecies))
		for i = 1:length(rMF.uniquespecies)
			y1 = X[c, i, :]
			y2 = Xe[c, i, :]
			pl = Gadfly.layer(x=rMF.tensorperiod.+1, y=y1, Gadfly.Geom.line(), Gadfly.Theme(line_width=2Gadfly.pt, default_color=parse(Colors.Colorant, scolors[i])))
			pd = Gadfly.layer(x=rMF.tensorperiod.+1, y=y2, Gadfly.Geom.line(), Gadfly.Theme(line_width=3Gadfly.pt, line_style=:dashdotdot, default_color=parse(Colors.Colorant, scolors[i])))
			pp = Gadfly.plot(pl..., pd..., Gadfly.Coord.Cartesian(xmin=2005, xmax=2017), Gadfly.Guide.xlabel(""), Gadfly.Guide.ylabel(rMF.uniquespecies[i]), Gadfly.Guide.title(w))
			Gadfly.draw(Gadfly.PNG("$figuredir/$prefix$w-$(rMF.uniquespecies[i]).png", 6Gadfly.inch, 3Gadfly.inch, dpi=300), pp)
			p[i] = Gadfly.plot(pl..., pd..., Gadfly.Coord.Cartesian(xmin=2005, xmax=2017), Gadfly.Guide.xlabel(""), Gadfly.Guide.ylabel(rMF.uniquespecies[i]), Gadfly.Theme(major_label_font_size=24Gadfly.pt, minor_label_font_size=18Gadfly.pt))
		end
		f = Compose.compose(Compose.context(0, 0, 1Compose.w, 0.001Compose.h),
			(Compose.context(), Compose.fill("gray"), Compose.fontsize(10Compose.pt), Compose.text(0.5, 1, w, Compose.hcenter, Compose.vtop)))
		pl = Gadfly.vstack(f, p...)
		Gadfly.draw(Gadfly.PNG("$figuredir/$prefix$w-allconc.png", 6Gadfly.inch, 3Gadfly.inch * length(rMF.uniquespecies), dpi=300), pl)
		Gadfly.draw(Gadfly.PDF("$figuredir/$prefix$w-allconc.pdf", 6Gadfly.inch, 3Gadfly.inch * length(rMF.uniquespecies)), pl)
	end

	info("Plot normalized species concentrations at each well (normalized by set) ...")
	pl = Vector{Any}(length(rMF.uniquespecies))
	pd = Vector{Any}(length(rMF.uniquespecies))
	for (c, w) in enumerate(rMF.wellnameorder)
		for i = 1:length(rMF.uniquespecies)
			Xn[c, i, :] ./= smax[i]
			y1 = X[c, i, :] / smax[i]
			y2 = Xe[c, i, :] / smax[i]
			pl[i] = Gadfly.layer(x=rMF.tensorperiod.+1, y=y1, Gadfly.Geom.line(), Gadfly.Theme(line_width=2Gadfly.pt, default_color=parse(Colors.Colorant, scolors[i])))
			pd[i] = Gadfly.layer(x=rMF.tensorperiod.+1, y=y2, Gadfly.Geom.line(), Gadfly.Theme(line_width=2Gadfly.pt, line_style=:dash,default_color=parse(Colors.Colorant, scolors[i])))
		end
		f = Gadfly.plot(pl..., pd..., Gadfly.Coord.Cartesian(xmin=2005, xmax=2017, ymin=0, ymax=1), Gadfly.Guide.xlabel(""), Gadfly.Guide.ylabel("Normalized concentrations"), Gadfly.Guide.title(w), Gadfly.Guide.manual_color_key("Species", rMF.uniquespecies, scolors[1:length(rMF.uniquespecies)]))
		Gadfly.draw(Gadfly.PNG("$figuredir/$prefix$w-normallconc.png", 6Gadfly.inch, 3Gadfly.inch, dpi=300), f)
		display(f)
		println()
	end

	# info("Plot normalized species concentrations at each well (normalized by well) ...")
	# pl = Vector{Any}(length(rMF.uniquespecies))
	# pd = Vector{Any}(length(rMF.uniquespecies))
	# for (c, w) in enumerate(rMF.wellnameorder)
	# 	for i = 1:length(rMF.uniquespecies)
	# 		y1 = (X[c, i, :] - sminw[c, i]) / (smaxw[c, i] - sminw[c, i])
	# 		y2 = (Xe[c, i, :] - sminw[c, i]) / (smaxw[c, i] - sminw[c, i])
	# 		pl[i] = Gadfly.layer(x=rMF.tensorperiod.+1, y=y1, Gadfly.Geom.line(), Gadfly.Theme(line_width=2Gadfly.pt, default_color=parse(Colors.Colorant, scolors[i])))
	# 		pd[i] = Gadfly.layer(x=rMF.tensorperiod.+1, y=y2, Gadfly.Geom.line(), Gadfly.Theme(line_width=2Gadfly.pt, line_style=:dash,default_color=parse(Colors.Colorant, scolors[i])))
	# 	end
	# 	f = Gadfly.plot(pl..., pd..., Gadfly.Coord.Cartesian(xmin=2005, xmax=2017, ymin=0, ymax=1), Gadfly.Guide.xlabel(""), Gadfly.Guide.ylabel("Normalized concentrations"), Gadfly.Guide.title(w), Gadfly.Guide.manual_color_key("Species", rMF.uniquespecies, scolors[1:length(rMF.uniquespecies)]))
	# 	Gadfly.draw(Gadfly.PNG("$figuredir/$prefix$w-normwellconc.png", 6Gadfly.inch, 3Gadfly.inch, dpi=300), f)
	# 	display(f)
	# 	println()
	# end

	info("NTFk results for w$(rMF.wellssetid) s$(rMF.speciessetid) y$(rMF.tensorperiod) $ns sources ...")
	pl = Vector{Any}(ns)
	for (c, w) in enumerate(rMF.wellnameorder)
		for i = 1:ns
			pl[i] = Gadfly.layer(x=rMF.tensorperiod.+1, y=W[c, i, :], Gadfly.Geom.line(), Gadfly.Theme(line_width=2Gadfly.pt, default_color=parse(Colors.Colorant, scolors[i])))
		end
		f = Gadfly.plot(pl..., Gadfly.Coord.Cartesian(xmin=2005, xmax=2017, ymin=0, ymax=1), Gadfly.Guide.xlabel(""), Gadfly.Guide.ylabel("Source mixing"), Gadfly.Guide.title(w), Gadfly.Guide.manual_color_key("Sources", sourcenames, scolors[1:ns]))
		Gadfly.draw(Gadfly.PNG("$figuredir/$prefix$w-mixing.png", 6Gadfly.inch, 3Gadfly.inch, dpi=300), f)
		display(f)
		println()
	end

	# irange = [1,2,5]
	# info("NTFk results for w$(rMF.wellssetid) s$(rMF.speciessetid) y$(rMF.tensorperiod) $ns sources; plot only sources $(irange) ...")
	# for (c, w) in enumerate(rMF.wellnameorder)
	# 	pl = Vector{Any}(0)
	# 	for i = irange
	# 		p = Gadfly.layer(x=rMF.tensorperiod.+1, y=W[c, i, :], Gadfly.Geom.line(), Gadfly.Theme(line_width=2Gadfly.pt, default_color=parse(Colors.Colorant, scolors[i])))
	# 		push!(pl, p)
	# 	end
	# 	f = Gadfly.plot(pl..., Gadfly.Coord.Cartesian(xmin=2005, xmax=2017, ymin=0), Gadfly.Guide.xlabel(""), Gadfly.Guide.ylabel("Source mixing"), Gadfly.Guide.title(w), Gadfly.Guide.manual_color_key("Sources", sourcenames[irange], scolors[irange]))
	# 	Gadfly.draw(Gadfly.PNG("$figuredir/$prefix$w-125mixing.png", 6Gadfly.inch, 3Gadfly.inch), f)
	# 	display(f)
	# 	println()
	# end
end

function getalltensorresults(X::Array, nsrange; retries=1000, resultdir="tensor-results", keyword="scale-", prefix="", quiet=true, isn=isnan.(X))
	Xbest = nothing
	Wbest = nothing
	Hbest = nothing
	for ns = nsrange
		if prefix == ""
			if isdefined(rMF, :wellssetid)
				filename = "$resultdir/$(rMF.case)-w$(rMF.wellssetid)-s$(rMF.speciessetid)-y$(rMF.tensoryearstep)-$keyword$ns-$retries-all.jld"
			else
				filename = "$resultdir/$(rMF.case)-$keyword$ns-$retries-all.jld"
			end
		else
			filename = "$resultdir/$(prefix)-$keyword$ns-$retries-all.jld"
		end
		if isfile(filename)
			fit, Wbest, Hbest, W, H = JLD.load(filename, "fit", "Wbest", "Hbest", "W", "H")
			NMFk.fixmixers!(X, Wbest)
			Xbest = NMFk.mixmatchcompute(X, Wbest, Hbest, isn)
			ofc = vecnorm(X[.!isn] .- Xbest[.!isn])
			println("$filename: $ns: fit $(minimum(fit)) ($ofc)")
			!quiet && for i = 1:length(W)
				NMFk.fixmixers!(X, W[i])
				isnw = isnan(W[i])
				Xe = NMFk.mixmatchcompute(X, W[i], H[i], isn)
				ofc = vecnorm(X[.!isn] .- Xe[.!isn])
				@show ofc, maximum(W[i][.!isnw]), maximum(H[i])
			end
		else
			warn("$filename is missing")
		end
	end
	return Xbest, Wbest, Hbest
end

function gettensorresults(X::Array, nsrange; retries=1000, resultdir="tensor-results", keyword="scale-", prefix="", quiet=true, isn=isnan.(X))
	Xe = nothing
	W = nothing
	H = nothing
	for ns = nsrange
		if prefix == ""
			if isdefined(rMF, :wellssetid)
				filename = "$resultdir/$(rMF.case)-w$(rMF.wellssetid)-s$(rMF.speciessetid)-y$(rMF.tensoryearstep)-$keyword$ns-$retries.jld"
			else
				filename = "$resultdir/$(rMF.case)-$keyword$ns-$retries.jld"
			end
		else
			filename = "$resultdir/$(prefix)-$keyword$ns-$retries.jld"
		end
		if isfile(filename)
			of, sil, aic, W, H = JLD.load(filename, "fit", "robustness", "aic", "W", "H")
			NMFk.fixmixers!(X, W)
			Xe = NMFk.mixmatchcompute(X, W, H, isn)
			ofc = vecnorm(X[.!isn] .- Xe[.!isn])
			println("$filename: $ns: fit $of ($ofc) silhouette $sil aic $aic")
		else
			warn("$filename is missing")
		end
	end
	return Xe, W, H
end

function getoldtensorresults(X::Array, nsrange; retries=1, resultdir="tensor-results", keyword="", prefix="", quiet=true, isn=isnan.(X))
	Xe = nothing
	W = nothing
	H = nothing
	for ns = nsrange
		if prefix == ""
			if isdefined(rMF, :wellssetid)
				filename = "$resultdir/$(rMF.case)-w$(rMF.wellssetid)-s$(rMF.speciessetid)-y$(rMF.tensoryearstep)-$keyword$ns-$retries.jld"
			else
				filename = "$resultdir/$(rMF.case)-$keyword$ns-$retries.jld"
			end
		else
			filename = "$resultdir/$(prefix)-$keyword$ns-$retries.jld"
		end
		if isfile(filename)
			of, W, H = JLD.load(filename, "OF", "W", "H")
			NMFk.fixmixers!(X, W)
			Xe = NMFk.mixmatchcompute(X, W, H, isn)
			ofc = vecnorm(X[.!isn] .- Xe[.!isn])
			println("$filename: $ns: fit $of ($ofc)")
		else
			warn("$filename is missing")
		end
	end
	return Xe, W, H
end

function roundsf(d::Number)
	f = floor(log10(d))
	rbase = convert(Int64, ceil(10^f/2))
	if f < 0
		rbase = 10
		sigdig = convert(Int64, -f)
	else
		sigdig = -1
	end
	round(d, sigdig, rbase)
end