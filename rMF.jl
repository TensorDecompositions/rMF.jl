module rMF

import NMFk
import JLD
import StatsBase
import DataStructures
import HypothesisTests
import Distributions
import SpatialAnalysis
import PyCall
import Gadfly
import Compose
import Colors
import PyPlot
import Images
import Compat
import Compat.AbstractString

@PyCall.pyimport matplotlib.patheffects as PathEffects

if isfile(Pkg.dir("Mads") * "/scripts/madsdisplay.jl")
	include(Pkg.dir("Mads") * "/scripts/madsdisplay.jl")
end

maxbuckets = 10
case = ""
casekeyword = ""
mixers = Array{Array{Float64, 2}}(maxbuckets)
buckets = Array{Array{Float64, 2}}(maxbuckets)
fitquality = Array{Float64}(maxbuckets)
robustness = Array{Float64}(maxbuckets)
aic = Array{Float64}(maxbuckets)
deltadependency = Array{Int64}(0)
dict_species = Dict()
uniquewells = Array{String}(0)
uniquespecies = Array{String}(0)
truebucket = Array{Float64}(0)
truedeltas = Array{Float64}(0)
truemixer = Array{Float64}(0)

include("rMFtranspose.jl")
include("rMFloaddata.jl")
include("rMFexecute.jl")
include("rMFgetresults.jl")

#=
dict_species = DataStructures.OrderedDict{AbstractString,AbstractString}(
	"Acidity or Alkalinity of a solution"=>"pH",
	"Alkalinity-CO3+HCO3"=>"CO3",
	"Chromium"=>"Cr",
	"Chromium-53/52"=>"δ53Cr",
	"Bromide"=>"Br-",
	"Calcium"=>"Ca",
	"Chloride"=>"Cl-",
	"Chlorate"=>"ClO3",
	"Perchlorate"=>"ClO4",
	"CHLORINE-36"=>"36Cl",
	"Chlorine-36/Chlorine Ratio (e-15)"=>"r36Cl",
	"Tritium"=>"3H",
	"Deuterium Ratio"=>"δ2H",
	"Oxygen-18/Oxygen-16 Ratio"=>"δ18O",
	"Magnesium"=>"Mg",
	"Nitrate-Nitrite as Nitrogen"=>"NO3",
	"Nitrogen-15/Nitrogen-14 Ratio(NO3)"=>"δ15N",
	"Nitrogen-15/Nitrogen-14 Ratio"=>"δ15N",
	"Oxygen-18/Oxygen-16 Ratio from Nitrate"=>"δ18O-NO3",
	"Sulfate"=>"SO4",
	"Sulfur-34/Sulfur-32 Ratio (SO4)"=>"δ34S-SO4",
	"Oxygen-18/Oxygen-16 Ratio from SO4"=>"δ18O-SO4",
	"Potassium"=>"K",
	"Rhenium"=>"Re",
	"Sodium"=>"Na",
	"Iodate"=>"IO3",
	"Iodine-129"=>"129I",
	"Iodine-129/Iodine Ratio (e-15)"=>"r129I",
	"Fraction Modern Carbon (de-normalized)"=>"r14C",
	"Dioxane[1,4-]"=>"Dioxane",
	"Acetaminophen"=>"Acetam",
	"Caffeine"=>"Caffe",
	"Sulfamethoxazole"=>"Sulfame")
JLD.save("data/cr-species.jld", "species", dict_species)

deltamixtures = DataStructures.OrderedDict{Any,Any}(
	"δ53Cr"=>ASCIIString["Cr"],
	"δ2H"=>Any[],
	"δ18O"=>Any[],
	"δ15N"=>ASCIIString["NO3"],
	"δ18O-NO3"=>ASCIIString["NO3"],
	"δ34S-SO4"=>ASCIIString["SO4"],
	"δ18O-SO4"=>ASCIIString["SO4"])
JLD.save("data/cr-stable-isotope-mixtures.jld", "deltamixtures", deltamixtures)

isotoperatios = DataStructures.OrderedDict{Any,Any}(
	"r36Cl"=>ASCIIString["Cl-"],
	"rI129"=>ASCIIString["I"])
JLD.save("data/cr-stable-isotope-ratios.jld", "isotoperatios", isotoperatios)

deltastandards = DataStructures.OrderedDict{Any,Any}(
	"δ53Cr"=>0.11339,
	"δ15N"=>0.003676,
	"δ18O-NO3"=>0.0.0020052,
	"δ34S-SO4"=>0.045005,
	"δ18O-SO4"=>0.0.0020052)
JLD.save("data/cr-stable-isotope-standards.jld", "deltastandards", deltastandards)
=#

info("rMF (Robust Matrix Factorization)")
info("")
info("Use `rMF.loaddata()` to get data")
info("Use `rMF.execute(5, retries=50)` to compute the results for the 5 bucket case with 50 reruns")
info("Use `rMF.getresults(5, retries=50)` to get the results for the 5 bucket case.")
info("Use `rMF.getresults(5:8, retries=50; brief=true)` to see the results for bucket cases 5 to 8.")
info("")
info("Use `rMF.loaddata()` to load different data sets:")
info("rMF.loaddata(20151202) - original fingerprint data set")
info("rMF.loaddata(20160202) - new fingerprint data set with pz wells")
info("""rMF.loaddata(20160202; wellsset="01", speciesset="13")""")
info("""rMF.loaddata("rdx-20160721")""")
info("""rMF.loaddata("test56s4")""")
info("""rMF.loaddata("test", nw=6, nc=4, ns=3)""")
info("")
# info("Test rMF")
# include(joinpath(Pkg.dir("rMF"), "test", "runtests.jl"))
info("Have fun ...")

rmfdir = splitdir(splitdir(Base.source_path())[1])[1]
cd(joinpath(rmfdir, "AquiferMixing"))
# cd(joinpath(Pkg.dir("rMF"), "AquiferMixing"))

end
