using LinearAlgebra, Printf, DelimitedFiles, Random
using JuMP, Gurobi, CPLEX
using Graphs, GraphPlot, SimpleWeightedGraphs
using Distances
using Statistics, Distributions, MultivariateStats
using Cairo, Compose, Meshes


const Optimizer = "CPLEX"  # ARGS[1]
if Optimizer == "Gurobi"
    @info "Using Gurobi"
    ENV["JULIA_DEBUG"] = Main
    # avoids checking the license everytime a model is created
    using Gurobi
    const GUROBI_ENV = Gurobi.Env()
else
    @info "Using CPLEX"
    using CPLEX
end
const ϵ = 0.001;
const TIMELIMIT = 3600;
const THREADS = 4;

include("data.jl")
include("generate_data.jl")
include("algos.jl")
include("X_and_cost.jl")
include("run.jl")