using LinearAlgebra, Printf, DelimitedFiles, Random
using JuMP, Gurobi, CPLEX
using Graphs, GraphPlot, SimpleWeightedGraphs
using Distances
using Statistics, Distributions, MultivariateStats
using Cairo, Compose, Meshes

#ENV["JULIA_DEBUG"] = Main

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
const THREADS = 4;
const output_flag = 0

function create_model()
    model = Model()
    if Optimizer == "Gurobi"
       set_optimizer(model, () -> Gurobi.Optimizer(GUROBI_ENV))
       set_optimizer_attribute(model, "OutputFlag", output_flag)
       set_optimizer_attribute(model, "TuneOutput", 1)
       set_optimizer_attribute(model, "TimeLimit", TIMELIMIT)
       set_optimizer_attribute(model, "NodefileStart", 0.5)
    else
       set_optimizer(model, () -> CPLEX.Optimizer())    
       set_optimizer_attribute(model, "CPX_PARAM_SCRIND", output_flag)
       set_optimizer_attribute(model,"CPX_PARAM_MIPDISPLAY", 0)
       set_optimizer_attribute(model, "CPX_PARAM_TILIM", TIMELIMIT)
       set_optimizer_attribute(model, "CPX_PARAM_NODEFILEIND", 2)
       set_optimizer_attribute(model, "CPX_PARAM_THREADS", THREADS)
       set_optimizer_attribute(model, "CPXPARAM_ScreenOutput", output_flag)
    end
    return model
 end

include("data.jl")
include("generate_data.jl")
include("algos.jl")
include("X_and_cost.jl")
include("run.jl")