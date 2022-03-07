abstract type Data end

#-------------------------------------------------------------------------------

struct Data_STP <: Data
  instance::String
  n::Int #number of nodes
  m::Int #number of edges
  g::SimpleGraph  # Underlying graph
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
  Δ::Float64
  E::Vector{Tuple{Int64,Int64}} #List the edges of the ground graph
  cost::Array{Array{Float64,2},2} # #cost[i,j][u_i,u_j] = distance between u_i ∈ U_i and u_j ∈ U_j
  c_center::Dict{Tuple{Int64,Int64},Float64}  # nominal distances, e.g. d(i,j)
  c_max::Dict{Tuple{Int64,Int64},Float64}  # maximum pairwise distances
  c_avg::Dict{Tuple{Int64,Int64},Float64}  # average pairwise distances

  function Data_STP(instance, n, m, g, from, to, δ⁻, δ⁺, t, t′, T, b, pos, U, nU, Δ, E, c_center)
    # cost in Steiner trees is given by Euclidean distances
    cost = Array{Array{Float64,2},2}(undef, n, n);
    for i in 1:n
      for j in i:n
        cost[i,j] = [norm(U[j][l]-U[i][k]) for l in 1:nU, k in 1:nU] 
        cost[j,i] = cost[i,j];
      end
    end

    # compute maximum and average pairwise costs
    c_max = Dict()
    for e in E c_max[e] = maximum(cost[e[1],e[2]]) end
    c_avg = Dict()
    for e in E c_avg[e] = mean(cost[e[1],e[2]]) end

    return new(instance, n, m, g, from, to, δ⁻, δ⁺, t, t′, T, b, pos, U, nU, Δ, E, cost, c_center, c_max, c_avg);
  end
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
  E::Vector{Tuple{Int64,Int64}} #List the edges of the ground graph
  cost::Array{Array{Float64,2},2} # #cost[i,j][u_i,u_j] = distance between u_i ∈ U_i and u_j ∈ U_j
  c_center::Dict{Tuple{Int64,Int64},Float64}  # nominal distances, e.g. d(i,j)
  c_max::Dict{Tuple{Int64,Int64},Float64}  # maximum pairwise distances
  c_avg::Dict{Tuple{Int64,Int64},Float64}  # average pairwise distances

  function Data_UFLP(instance, n, m, g, I, J, nU, U, p, E, c_center) 
    # cost in UFLP is given by Euclidean distances
    cost = Array{Array{Float64,2},2}(undef, n, n);
    for i in 1:n
      for j in i:n
        cost[i,j] = [norm(U[j][l]-U[i][k]) for l in 1:nU, k in 1:nU] 
        cost[j,i] = cost[i,j];
      end
    end

    # compute maximum and average pairwise costs
    c_max = Dict()
    for e in E c_max[e] = maximum(cost[e[1],e[2]]) end
    c_avg = Dict()
    for e in E c_avg[e] = mean(cost[e[1],e[2]]) end

    return new(instance, n, m, g, I, J, nU, U, p, E, cost, c_center, c_max, c_avg)
  end
end

"""
  Data_clustering
"""
struct Data_clustering <: Data
  instance::String  # name of the dataset
  n::Int  # number of vertices
  dimension::Int  # dimension of the metric space
  nU::Vector{Int64}  # number of points in each uncertainty set
  U::Vector{Vector{Vector{Float64}}}  # coordinates of the points in the uncertainty sets of each node
  K::Int  # number of clusters targeted in a balanced clustering
  E::Vector{Tuple{Int64,Int64}}  # List the edges of the ground graph
  cost::Array{Array{Float64,2},2} # pairwise squared squared euclidean distances between each pair of uncertainty sets
  c_center::Dict{Tuple{Int64,Int64},Float64}  # nominal distances, e.g. d(i,j)
  c_max::Dict{Tuple{Int64,Int64},Float64}  # maximum pairwise distances
  c_avg::Dict{Tuple{Int64,Int64},Float64}  # average pairwise distances

  function Data_clustering(instance, n, dim, nU, U, K)
    E = [(i,j) for i in 1:n for j in i+1:n]

    # cost in clustering is given by squared Euclidean distances
    cost = Array{Array{Float64,2},2}(undef, n, n);
    for i in 1:n
      for j in i:n
        cost[i,j] = zeros(nU[i],nU[j])
        for k in 1:nU[i]
          for l in 1:nU[j]
            cost[i,j][k,l] = sum((U[j][l]-U[i][k]).^2)
          end
        end
        cost[j,i] = cost[i,j];
      end
    end

    # compute the paiwise squared distances between the barycenters of U
    barycenter = Vector{Vector{Float64}}();
    for i ∈ 1:n
        push!(barycenter, sum(U[i][k] for k in 1:nU[i])./nU[i])
    end
    c_center = Dict()
    for e in E
      c_center[e] = sum((barycenter[e[2]]-barycenter[e[1]]).^2)
    end

    # compute maximum and average pairwise costs
    c_max = Dict()
    for e in E c_max[e] = maximum(cost[e[1],e[2]]) end
    c_avg = Dict()
    for e in E c_avg[e] = mean(cost[e[1],e[2]]) end

    return new(instance, n, dim, nU, U, K, E, cost, c_center, c_max, c_avg)
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

