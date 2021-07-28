using JuMP, Gurobi, CPLEX
using LightGraphs, GraphPlot, SimpleWeightedGraphs
using LinearAlgebra
using Printf, DelimitedFiles
using Random, MultivariateStats

const Optimizer = "CPLEX"  # ARGS[1]
if Optimizer == "Gurobi"
    @info "Using Gurobi"
    # avoids checking the license everytime a model is created
    using Gurobi
    const GUROBI_ENV = Gurobi.Env()
else
    @info "Using CPLEX"
    using CPLEX
end
const ϵ = 0.001;
const TIMELIMIT = 7200;
const THREADS = 8;

include("data.jl")
include("algos.jl")
include("X_and_cost.jl")


#------

function printres(time, res, seed, algo, data::Data_STP)
    f = open("res/Steiner/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.Δ) $(data.nU) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end

#------

function printres(time, res, seed, algo, data::Data_UFLP)
    f = open("res/UFLP/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.n) $(data.m) $(data.nU) $(length(data.I)) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end

#------

function run_steiner()
    Random.seed!(1)
	seed = 0

    @warn "warming up ..."
    data = read_data_STP("data/Steiner/small.stp", 0, 1)
    @info "CP for $(data.instance) with Δ=$(data.Δ) and $(data.nU) extreme points"
    exact(data)
    heuristic_det(data)
    heuristic_dmax(data)
    heuristic_adr(data)

    @warn "now looping over all problems"

	"NOTE: Instances small and small_i correspond to the instances format_i from the paper"
    #folder = "small_i/"
	folder = "P6E/"
    #INSTANCES = readdir("data/Steiner/"*folder)
    INSTANCES =["p619.stp"] #,"p620.stp","p621.stp"]
    for instance in INSTANCES, nU in [10], Δ in 0.1
    #for nU in [5,10,20], Δ in [0.1,0.5,1], seed in 1:20, size in 1:2
        #@warn "seed number $seed"
        data = read_data_STP("data/Steiner/"*folder*instance,Δ,nU)
        #data = create_small_STP(size,Δ,nU)

		time_exact = @elapsed res_exact = exact(data);
	    time_worst = @elapsed res_worst = heuristic_dmax(data);
	    time_center = @elapsed res_center = heuristic_det(data);
		#time_adr = @elapsed res_adr = heuristic_adr(data);

		printres(time_exact, res_exact, seed, folder*"exact", data);
	    printres(time_worst, res_worst, seed, folder*"worst", data);
	    printres(time_center, res_center, seed, folder*"center", data);
		#printres(time_adr, res_adr, seed, folder*"adr", data);

        flush(stdout)
    end
end

#------

function run_UFLP()
	@info "Warming up ..."
	n, m, cardI, nU, p, seed = 10, 20, 1, 1, 1, 1
	data = build_UFLP(n, m, cardI, nU, p, seed)
	res_exact = exact(data);
	res_worst = heuristic_dmax(data);
	res_center = heuristic_det(data);
	@info "... finished!"

	n, cardI, p = 80, 15, 3
	maxU = Int64(floor(n/cardI))
	for seed in 1:20, m in [n, round(n*1.5), round(n*2)], nU in [3,4,5]
		@warn "($n, $m, $cardI, $nU, $seed, $p)"
		data = build_UFLP(n, m, cardI, nU, seed, p)

		time_exact = @elapsed res_exact = exact(data);
	    time_worst = @elapsed res_worst = heuristic_dmax(data);
	    time_center = @elapsed res_center = heuristic_det(data);

		printres(time_exact, res_exact, seed, "exact", data);
	    printres(time_worst, res_worst, seed, "worst", data);
	    printres(time_center, res_center, seed, "center", data);
	end
end

#------

#run_UFLP()
run_steiner()
