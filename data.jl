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

#--------
# Generation of instances along the lines of M. Daskin. Genrand2: A random network generator.
# Department of Industrial Engineering and Management Sciences, Northwestern University,
# Evanston, IL, USA, 1993
function build_UFLP(n, m, cardI, nU, seed, p)
	Random.seed!(seed)
	positions = rand(n,2)
	# 1- ensures connexity by adding a min-cost spanning tree, favorize short edges to emulate transportation networks
	sources = [i for i in 1:n for j in (i+1):n]
	destinations = [j for i in 1:n for j in (i+1):n]
	g_MST = complete_graph(n)
	dst = [norm(positions[i,:] - positions[j,:])^2 for i in 1:n, j in 1:n]
	MST = kruskal_mst(g_MST, dst)
	g = SimpleWeightedGraph(n)
	for e in MST add_edge!(g,e.src,e.dst,norm(positions[e.src,:] - positions[e.dst,:])) end

	#draw(PDF("UFLP0.pdf", 16cm, 16cm), gplot(g, positions[:,1], positions[:,2], nodelabel=1:nv(g)))

	# 2- add further edges not crossing existing ones to emulate transportation networks
	dst = [1/norm(positions[i,:] - positions[j,:])^2 for i in 1:n for j in (i+1):n]
	edges = [(i,j) for i in 1:n for j in (i+1):n]
	E_MST = [(e.src,e.dst) for e in MST]
	indexes_not_in_MST = findall([!in(e,E_MST) for e in edges])
	edges = edges[indexes_not_in_MST]
	dst = dst[indexes_not_in_MST]
	# removing edges crossing MST
	indexes_crossing_MST = []
	for i in 1:length(edges)
		e = edges[i]
		s = Segment((positions[e[1],1],positions[e[1],2]), (positions[e[2],1],positions[e[2],2]))
		for j in 1:length(E_MST)
			e′ = E_MST[j]
			s′ = Segment((positions[e′[1],1],positions[e′[1],2]), (positions[e′[2],1],positions[e′[2],2]))
			if s ∩ s′ != nothing
				# two edges sharing an extremeity are not considered crossing
				point = coordinates(s ∩ s′)
				dst_to_endpoints = min(norm(point-positions[e[1],:]),norm(point-positions[e[2],:]),norm(point-positions[e′[1],:]),norm(point-positions[e′[2],:]))
				if dst_to_endpoints > 0.0001
					push!(indexes_crossing_MST,i)
					break
				end
			end
		end
	end
	edges = edges[setdiff(1:length(edges),indexes_crossing_MST)]
	dst = dst[setdiff(1:length(dst),indexes_crossing_MST)]
	while ne(g) < m
		if isempty(edges)
			@error "Not possible to generate the graph"
			exit()
		end
		M = length(dst)
		# normalize the probabilities
		probas = dst/sum(dst)
		# define the cumulative distribution
		cumulatives = zeros(M+1)
		cumulatives[2:end] = [sum(probas[1:i]) for i in 1:M]
		cumulatives[end] = 1
		new_rand = rand()
		for i in 1:M
			if cumulatives[i] ≤ new_rand && new_rand ≤ cumulatives[i+1]
				new_edge = i
				e′ = edges[i]
				add_edge!(g,e′[1],e′[2],norm(positions[e′[1],:] - positions[e′[2],:]))
				edges = edges[setdiff(1:M,i)]
				dst = dst[setdiff(1:M,i)]
				# removing edges crossing e′
				indexes_crossing_e′ = []
				for j in 1:length(edges)
					e = edges[j]
					s = Segment((positions[e[1],1],positions[e[1],2]), (positions[e[2],1],positions[e[2],2]))
					s′ = Segment((positions[e′[1],1],positions[e′[1],2]), (positions[e′[2],1],positions[e′[2],2]))
					if s ∩ s′ != nothing
						# two edges sharing an extremeity are not considered crossing
						point = coordinates(s ∩ s′)
						dst_to_endpoints = min(norm(point-positions[e[1],:]),norm(point-positions[e[2],:]),norm(point-positions[e′[1],:]),norm(point-positions[e′[2],:]))
						if dst_to_endpoints > 0.0001
							push!(indexes_crossing_e′,j)
						end
					end
				end
				edges = edges[setdiff(1:end,indexes_crossing_e′)]
				dst = dst[setdiff(1:end,indexes_crossing_e′)]
				break
			end
		end
	end
	instance = "random_$(n)_$(m)_$(nU)"
	U = Vector{Vector{Int64}}()
	dst = [norm(positions[i,:] - positions[j,:]) for i in 1:n, j in 1:n]
	mean_dst = sum(dst)/(n*(n-1)/2)
	dst = Matrix{Float64}(undef,n,n)
	for i in 1:n
		ds = dijkstra_shortest_paths(g, i)
		neighbours = collect(1:n)[sortperm(ds.dists)]
		neighbours = neighbours[1:nU] # take the nU clostest neigbours
		push!(U,neighbours)
		dst[i,:] = ds.dists
	end
	V = collect(1:n)
	I = []
	while length(I) < cardI
		if isempty(V)
			@error "Not possible to generate the full set I"
			exit()
		end
		i = rand(V)
		push!(I,i)
		V = setdiff(V,U[i])
	end
	I = sort(I)
	J = setdiff(1:n,I)
	d = fill(-1.0,n,maximum(length.(U)),n,maximum(length.(U)))
	for i in 1:n, u_i in 1:length(U[i]), j in 1:n, u_j in 1:length(U[j])
		d[i,u_i,j,u_j] = dst[U[i][u_i],U[j][u_j]]
	end
	E = Vector{Tuple{Int64,Int64}}()
	d⁰ = Dict()
	for i in I, j in J
		push!(E,(i,j))
		d⁰[(i,j)] = dst[i,j]
	end
	#draw(PDF("UFLP.pdf", 16cm, 16cm), gplot(g, positions[:,1], positions[:,2], nodelabel=1:nv(g)))
	#@info I
	return Data_UFLP(instance, n, m, g, I, J, nU, U, p, d, E, d⁰)
end

#-------------------------------------------------------------------------------

function read_data_STP(instance,Δ,nU)
  datafile = readdlm(instance)
  line = findfirst(datafile.=="Nodes")[1]
  n = datafile[line, 2]
  line = findfirst(datafile.=="Edges")[1]
  m = datafile[line, 2]
  from = []
  to = []
  δ⁺=Vector{Vector{Int64}}()
  δ⁻=Vector{Vector{Int64}}()
  E = Vector{Tuple{Int64,Int64}}()

  for i in 1:n
    push!(δ⁺,[])
    push!(δ⁻,[])
  end
  for e in 1:m
    push!(from,datafile[line+e,2])
    push!(to,datafile[line+e,3])
    push!(δ⁺[from[e]],e)
    push!(δ⁻[to[e]],e)
    push!(δ⁺[to[e]],e+m)
    push!(δ⁻[from[e]],e+m)
	push!(E,(from[e],to[e]))
  end
  line = findfirst(datafile.=="Terminals")[1]
  t = datafile[line, 2]
  t′ = t-1
  T = []
  for i in 1:t
    push!(T,datafile[line+i,2])
  end
  b = zeros(n,t′)
  for tt in 1:t′
    b[T[tt],tt] = -1
    b[T[t],tt] = 1
  end
  line = findfirst(datafile.=="Coordinates")
  pos = Vector{Vector{Int64}}()
  if line != nothing
    line = line[1]
    for i in 1:n
      push!(pos,[datafile[line+i,3],datafile[line+i,4]])
    end
  else
    @info "Positions not included in data file => computed through a simple variant of MDS-MAP"
    line = findfirst(datafile.=="Edges")[1]
    weights = []
    for e in 1:m push!(weights,datafile[line+e,4]) end
    g = SimpleWeightedGraph(convert(Array{Int64,1},from), convert(Array{Int64,1},to), convert(Array{Float64,1},weights))
    D = zeros(n,n)
    for i in 1:n
      D[i,:] = dijkstra_shortest_paths(g, i).dists
    end
    positions = round.(transform(fit(MDS, D, maxoutdim=2, distances=true)))
    for i in 1:n
      push!(pos,[positions[1,i],positions[2,i]])
    end
  end
  distances = [norm(pos[i]-pos[j]) for i in 1:n, j in 1:n]
  distance_mean = sum(distances)/(n*(n-1))
  U = Vector{Vector{Vector{Float64}}}()
  for i in 1:n
    radius = rand() * Δ * distance_mean
    push!(U, Vector{Vector{Float64}}())
    for k in 1:nU
      push!(U[i],round.([pos[i][1]+radius*cos(2π*k/nU),pos[i][2]+radius*sin(2π*k/nU)]))
    end
  end
  d = zeros(n,nU,n,nU)
  for i in 1:n, k in 1:nU, j in 1:n, l in 1:nU
    d[i,k,j,l] = norm(U[i][k]-U[j][l])
  end
  d⁰ = Dict()
  for e in E
	d⁰[e] = distances[e[1],e[2]]
  end
  g = SimpleGraph(n)
  for e in E add_edge!(g, e[1], e[2]) end
  return Data_STP(instance,n,m,g,from,to,δ⁻,δ⁺,t,t′,T,b,pos,U,nU,d,Δ,E,d⁰)
end

#-------------------------------------------------------------------------------

function create_small_STP(dim,Δ,nU)
  instance = "small_$dim"
  n = 2+5*dim
  m = 1+8*dim
  T = []
  t = 1+2*dim
  t′ = t-1
  pos = Vector{Vector{Int64}}()
  b = zeros(n,t′)
  from = []
  to = []
  δ⁺=Vector{Vector{Int64}}()
  δ⁻=Vector{Vector{Int64}}()
  E = Vector{Tuple{Int64,Int64}}()
  for i in 1:n
    push!(δ⁺,[])
    push!(δ⁻,[])
  end

  push!(T,1)
  push!(from,1)
  push!(to,2)
  push!(pos,[55, 5])
  push!(pos,[105, 5])
  for layer in 1:dim
    term = 5*(layer-1)
    push!(T,term+5)
    push!(T,term+6)
    push!(from,term+1)
    push!(to,term+3)
    push!(from,term+2)
    push!(to,term+4)
    push!(from,term+3)
    push!(to,term+4)
    push!(from,term+2)
    push!(to,term+5)
    push!(from,term+3)
    push!(to,term+6)
    push!(from,term+4)
    push!(to,term+7)
    push!(from,term+5)
    push!(to,term+7)
    push!(from,term+6)
    push!(to,term+7)
    yterm = 90*(layer-1)
    push!(pos,[30,yterm+50])
    push!(pos,[80,yterm+50])
    push!(pos,[130,yterm+50])
    push!(pos,[55,yterm+95])
    push!(pos,[105,yterm+95])
  end

  for e in 1:m
    push!(δ⁺[from[e]],e)
    push!(δ⁻[to[e]],e)
    push!(δ⁺[to[e]],e+m)
    push!(δ⁻[from[e]],e+m)
	push!(E,(from[e],to[e]))
  end
  for tt in 1:t′
    b[T[tt],tt] = -1
    b[T[t],tt] = 1
  end
  distances = [norm(pos[i]-pos[j]) for i in 1:n, j in 1:n]
  distance_mean = sum(distances)/(n*(n-1))
  U = Vector{Vector{Vector{Float64}}}()
  for i in 1:n
    radius = rand() * Δ * distance_mean
    push!(U, Vector{Vector{Float64}}())
    for k in 1:nU
      push!(U[i],round.([pos[i][1]+radius*cos(2π*k/nU),pos[i][2]+radius*sin(2π*k/nU)]))
    end
  end
  d = zeros(n,nU,n,nU)
  for i in 1:n, k in 1:nU, j in 1:n, l in 1:nU
    d[i,k,j,l] = norm(U[i][k]-U[j][l])
  end
  d⁰ = Dict()
  for e in E
	d⁰[e] = distances[e[1],e[2]]
  end
  g = SimpleGraph(n)
  for e in E add_edge!(g, e[1], e[2]) end
  return Data_STP(instance,n,m,g,from,to,δ⁻,δ⁺,t,t′,T,b,pos,U,nU,d,Δ,E,d⁰)
end
