using Distributions, LinearAlgebra
using Distances
using LightGraphs, GraphPlot, Cairo, Compose, Meshes

abstract type Data end

#-------------------------------------------------------------------------------

struct Data_STP <: Data
  instance::String
  n::Int #number of nodes
  m::Int #number of edges
  g::SimpleGraph #Underlying graph
  from::Vector{Int64}
  to::Vector{Int64}
  δ⁻::Vector{Vector{Int64}}
  δ⁺::Vector{Vector{Int64}}
  t::Int
  t′::Int # =t-1
  T::Vector{Int64}
  b::Matrix{Int64}
  pos::Vector{Vector{Int64}} # nominal coordinates of each node
  U::Vector{Vector{Vector{Float64}}} # coordinates of the points in the uncertainty sets of each node
  nU::Int64 # number of points in each uncertainty set
  d::Array{Float64,4} #d(i,u_i,j,u_j) = distance between u_i ∈ U_i and u_j ∈ U_j
  Δ::Float64

  E::Vector{Tuple{Int64,Int64}} #List the edges of the ground graph
  d⁰::Dict #nominal distances, e.g. d(i,j)

  Data_STP(instance, n, m, g, from, to, δ⁻, δ⁺, t, t′, T, b, pos, U, nU, d, Δ, E, d⁰) =
  new(instance, n, m, g, from, to, δ⁻, δ⁺, t, t′, T, b, pos, U, nU, d, Δ, E, d⁰);
end

#-------------------------------------------------------------------------------

struct Data_UFLP <: Data
  instance::String #name of the instance
  n::Int #number of nodes
  m::Int #number of edges
  g::SimpleWeightedGraph #Underlying graph
  I::Vector{Int64} #indexes of clients
  J::Vector{Int64} #indexes of facilities
  nU::Int64 #size of U[i]
  U::Vector{Vector{Int64}} #indexes of nodes in U_i for each i in V
  p::Int #number of facilities to be installed
  d::Array{Float64,4} #d(i,u_i,j,u_j) = distance between U_i[u_i] and U_j[u_j]
  E::Vector{Tuple{Int64,Int64}} #List the edges of the ground graph
  d⁰::Dict #nominal distances, e.g. d(i,j)

  Data_UFLP(instance, n, m, g, I, J, nU, U, p, d, E, d⁰) = new(instance, n, m, g, I, J, nU, U, p, d, E, d⁰)
end

"""
  Data_clustering
"""
struct Data_clustering <: Data
  instance::String  # name of the dataset
  n::Int  # number of vertices
  dimension::Int  # dimension of the metric space
  nU::Vector{Int64}  # number of points in each uncertainty set
  U::Vector{Array{Float64,2}}  # uncertainty set of each vertex
  K::Int  # number of clusters targeted in a balanced clustering
  E::Vector{Pair{Int64}}  # List the edges of the ground graph
  d²::Array{Array{Float64,2},2} # paiwise squared squared euclidean distances between each pair of uncertainty sets
  d⁰::Dict  # nominal distances, e.g. d(i,j)

  function Data_clustering(instance, n, d, nU, U, K, d²)
    E = [(i,j) for i in 1:n, j in i+1:n]

    # compute the paiwise squared distances between the barycenters of U
    barycenter = Array{Float64,2}(undef,2,0);
    for i in V
       barycenter = [barycenter sum(data.U[i], dims=2)./data.nU[i]]
    end
    d⁰ = Dict()
    for e in E
      d⁰[e] = sum((barycenter[:,j]-barycenter[:,i]).^2)
    end

    return new(instance, n, d, nU, U, K, E, d², d⁰)
  end
end

"""
  Data_interval
"""
struct Data_interval <: Data
  instance::String  # name of the dataset
  n::Int  # number of vertices
  dimension::Int  # dimension of the metric space, i.e., number of attributes
  lb::Array{Float64, 2}  # lower bound of each attribute of each vertex
  ub::Array{Float64, 2}  # upper bound of each attribute of each vertex
  K::Int  # number of clusters targeted in a balanced clustering

  Data_interval(instance, n, dim, lb, ub, K) = new(instance, n, dim, lb, ub, K)
end

