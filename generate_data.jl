"""#--------
# Generation of instances along the lines of M. Daskin. Genrand2: A random network generator.
# Department of Industrial Engineering and Management Sciences, Northwestern University,"""
# Evanston, IL, USA, 1993
function build_UFLP(n, m, cardI, nU, seed, p)
  Random.seed!(seed)
  positions = rand(n,2)
  # 1- ensures connectivity by adding a min-cost spanning tree, favorize short edges to emulate transportation networks
  sources = [i for i in 1:n for j in (i+1):n]
  destinations = [j for i in 1:n for j in (i+1):n]
  g_MST = complete_graph(n)
  dst = [norm(positions[i,:] - positions[j,:])^2 for i in 1:n, j in 1:n]
  MST = kruskal_mst(g_MST, dst)
  g = SimpleWeightedGraph(n)
  cost = Dict()
  for e in MST 
    add_edge!(g,e.src,e.dst, norm(positions[e.src,:] - positions[e.dst,:]))
    cost[(e.src,e.dst)] = norm(positions[e.src,:] - positions[e.dst,:])
  end

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
      if s ∩ s′ !== nothing
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
          if s ∩ s′ !== nothing
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
  E = Vector{Tuple{Int64,Int64}}()
  c_center = Dict()
  for i in I, j in J
    push!(E,(i,j))
    c_center[(i,j)] = dst[i,j]
  end
  #draw(PDF("UFLP.pdf", 16cm, 16cm), gplot(g, positions[:,1], positions[:,2], nodelabel=1:nv(g)))
  #@info I
  return Data_UFLP(instance, n, m, g, I, J, nU, U, p, E, c_center)
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

  # TODO: this is really heavy for nothing, class Graph provides everything
  for i in 1:n
    push!(δ⁺,[])
    push!(δ⁻,[])
  end
  g = SimpleGraph(n)
  for e in 1:m
    push!(from,datafile[line+e,2])
    push!(to,datafile[line+e,3])
    push!(δ⁺[from[e]],e)
    push!(δ⁻[to[e]],e)
    push!(δ⁺[to[e]],e+m)
    push!(δ⁻[from[e]],e+m)
    push!(E,(from[e],to[e]))
    add_edge!(g, from[e], to[e])
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
  if line !== nothing
    line = line[1]
    for i in 1:n
      push!(pos,[datafile[line+i,3],datafile[line+i,4]])
    end
  else
    @info "Positions not included in data file => computed through a simple variant of MDS-MAP"
    line = findfirst(datafile.=="Edges")[1]
    weights = []
    for e in 1:m 
      push!(weights,datafile[line+e,4])
    end
    g_disjkstra = SimpleWeightedGraph(convert(Array{Int64,1},from), convert(Array{Int64,1},to), convert(Array{Float64,1},weights))
    D = zeros(n,n)
    for i in 1:n
      D[i,:] = dijkstra_shortest_paths(g_disjkstra, i).dists
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

  c_center = Dict()
  for e in E
    c_center[e] = distances[e[1],e[2]]
  end
  return Data_STP(instance,n,m,g,from,to,δ⁻,δ⁺,t,t′,T,b,pos,U,nU,Δ,E,c_center)
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
  c_center = Dict()
  for e in E
  c_center[e] = distances[e[1],e[2]]
  end
  g = SimpleGraph(n)
  for e in E add_edge!(g, e[1], e[2]) end
  return Data_STP(instance,n,m,g,from,to,δ⁻,δ⁺,t,t′,T,b,pos,U,nU,Δ,E,c_center)
end

"""
  build_gaussian_clustering
  Generation of synthetic bivariate gaussian datasets reproducing the procedure described in De Carvalho and Lechevallier (2009) "Partitional Clustering  Algorithms for Symbolic Interval Data Based on Single Adaptive Distances."
"""
function build_gaussian_clustering(name::String, n_per_cluster, interval_length, μ1, μ2, σ1², σ2², ρ12)
  K = length(μ1);
  # build bivariate gaussian distributions
  G = Vector{MultivariateDistribution}();
  for i in 1:K
    μ = [μ1[i];μ2[i]];
    σ1 = √σ1²[i];
    σ2 = √σ1²[i];
    Σ = [σ1^2 ρ12[i]*σ1*σ2 ; ρ12[i]*σ1*σ2 σ2^2];
    push!(G, MvNormal(μ, Σ));
  end

  # draw the centers of the uncertainty sets from the gaussian distributions
  # and then build each set as a rectangle with dimensions drawn uniformly in
  # an input interval
  U = Vector{Vector{Vector{Float64}}}();
  for i in 1:n_per_cluster
    for k in 1:K
      push!(U, Vector{Vector{Float64}}());
      z = rand(G[k]);
      γ = [1 + interval_length*rand() ; 1 + interval_length*rand()];
      push!(U[end], [z[1] - γ[1]/2.0 ; z[2] - γ[2]/2.0]);
      push!(U[end], [z[1] - γ[1]/2.0 ; z[2] + γ[2]/2.0]);
      push!(U[end], [z[1] + γ[1]/2.0 ; z[2] - γ[2]/2.0]); 
      push!(U[end], [z[1] + γ[1]/2.0 ; z[2] + γ[2]/2.0]);
    end
  end
  n = K * n_per_cluster;
  nU = 4 * Int.(ones(n));

  return Data_clustering(name, K * n_per_cluster, 2, nU, U, K);
end

"""
  build_together_clustering
 
  Generation of synthetic bivariate gaussian datasets reproducing the procedure described in De Carvalho and Lechevallier (2009) "Partitional Clustering Algorithms for Symbolic Interval Data Based on Single Adaptive Distances."
"""
function build_together_clustering(n_per_cluster, interval_length, μ1, μ2)
  K = length(μ1);

  # draw the centers of the uncertainty sets from the gaussian distributions
  # and then build each set as a rectangle with dimensions drawn uniformly in
  # an input interval
  U = Vector{Vector{Vector{Float64}}}();
  for i in 1:n_per_cluster
    for k in 1:K
      push!(U, Vector{Vector{Float64}}());
      z = [μ1[k];μ2[k]];
      γ = [1 + interval_length*rand() ; 1 + interval_length*rand()];
      push!(U[end], [z[1] - γ[1]/2.0 ; z[2] - γ[2]/2.0]);
      push!(U[end], [z[1] - γ[1]/2.0 ; z[2] + γ[2]/2.0]);
      push!(U[end], [z[1] + γ[1]/2.0 ; z[2] - γ[2]/2.0]); 
      push!(U[end], [z[1] + γ[1]/2.0 ; z[2] + γ[2]/2.0]);
    end
  end
  n = K * n_per_cluster;
  nU = 4 * Int.(ones(n));

  return Data_clustering("synthetic_clustering", K * n_per_cluster, 2, nU, U,
    K);
end

#-------------------------------------------------------------------------------
function add_dim(box::Vector{Vector{Float64}}, dim::Int, xmin::Vector{Float64}, xmax::Vector{Float64})
  nU = length(box)
  for k ∈ 1:nU
    box[k][dim] = xmin[dim]
    push!(box, copy(box[k]))
    box[end][dim] = xmax[dim]
  end
end
function create_box(U0::Vector{Float64}, dims::Vector{Int}, xmin::Vector{Float64}, xmax::Vector{Float64})
  box = Vector{Vector{Float64}}()
  push!(box, U0)
  for dim ∈ dims
    add_dim(box, dim, xmin, xmax)
  end
  return box
end

function read_balanced_clustering(instance,p_missing,K)
  coords = readdlm("data/BalancedClustering/"*instance)
  coords = coord[1:50,:]  # won't be solved for more than 50 nodes
  (n,dim) = size(coords)

  # Build uncertainty sets
  ## simulate missing data by p_missing percents components in each dimension
  ## missing data is replaced with intervals with extremities equal to the minimum and maximum observed values
  nb_missing = floor(Int, p_missing/100.0 * n)
  xmin = vec(minimum(coords, dims=1));
  xmax = vec(maximum(coords, dims=1));

  ## first draw the missing fields of each vertex
  missing_dims = Vector{Vector{Int}}()
  for i ∈ 1:n
    push!(missing_dims, Vector{Int}())
  end
  for k ∈ 1:dim
    missing = Vector{Int}()
    while length(missing) < nb_missing
      missing = unique(rand(1:n, 2*nb_missing))
    end
    for i ∈ missing[1:nb_missing]
      push!(missing_dims[i], k)
    end
  end

  ## then, create the corresponding uncertainty sets
  nU = Vector{Int}()
  U = Vector{Vector{Vector{Float64}}}()
  for i ∈ 1:n
    push!(nU, 2^(length(missing_dims[i])))
    push!(U, create_box(coords[i,:],missing_dims[i],xmin,xmax))
  end
  return Data_clustering(instance, n, dim, nU, U, K);
end

function read_cars(K::Int = 4, dim::Int = 2)
  coords = readdlm("data/IntervalClustering/cars.txt")
  (n,nvals) = size(coords)
  dim_init = round(Int, nvals/2)

  xmin = zeros(n,dim_init)
  xmax = zeros(n,dim_init)
  centers = zeros(n,dim_init)
  for i ∈ 1:n
    for j ∈ 1:dim_init
      xmin[i,j] = coords[i,2*j-1]
      xmax[i,j] = coords[i,2*j]
      centers[i,j] = (xmin[i,j]+xmax[i,j])/2
    end
  end
  centers_std = (centers .- mean(centers, dims=1))./std(centers, dims=1)
  xmin_std = (xmin .- mean(centers, dims=1))./std(centers, dims=1)
  xmax_std = (xmax .- mean(centers, dims=1))./std(centers, dims=1)

  M = fit(PCA, transpose(centers_std); maxoutdim=dim, pratio=1.0)
  xmin_pca = Array{Float64}(transpose(transform(M, transpose(xmin_std))))
  xmax_pca = Array{Float64}(transpose(transform(M, transpose(xmax_std))))
  centers_pca = Array{Float64}(transpose(transform(M, transpose(centers_std))))

  ## then, create the corresponding uncertainty sets
  nU = Vector{Int}()
  U = Vector{Vector{Vector{Float64}}}()
  for i ∈ 1:n
    push!(nU, 2^dim)
    push!(U, create_box(centers_pca[i,:], collect(1:2),xmin_pca[i,:],xmax_pca[i,:]))
  end

  return U
end

function read_meteo_data(K::Int = 6, dim::Int = 2)
  tmin = readdlm("data/IntervalClustering/nrmmin.txt")
  tmax = readdlm("data/IntervalClustering/nrmmax.txt")
  tavg = readdlm("data/IntervalClustering/nrmavg.txt")
  pcp = readdlm("data/IntervalClustering/nrmpcp.txt")


  tavg_std = (tavg .- mean(tavg, dims=1))./std(tavg, dims=1)
  tmin_std = (tmin .- mean(tavg, dims=1))./std(tavg, dims=1)
  tmax_std = (tmax .- mean(tavg, dims=1))./std(tavg, dims=1)
  pcp_std = (pcp .- mean(pcp, dims=1))./std(pcp, dims=1)


  M = fit(PCA, transpose(tavg_std); maxoutdim=dim, pratio=1.0)
  tmin_pca = Array{Float64}(transpose(transform(M, transpose(tmin_std))))
  tmax_pca = Array{Float64}(transpose(transform(M, transpose(tmax_std))))
  tavg_pca = Array{Float64}(transpose(transform(M, transpose(tavg_std))))

  M = fit(PCA, transpose(pcp_std); maxoutdim=1, pratio=1.0)
  pcp_pca = Array{Float64}(transpose(transform(M, transpose(pcp_std))))

  ## then, create the corresponding uncertainty sets
  nU = Vector{Int}()
  n= size(tmin)[1]
  U = Vector{Vector{Vector{Float64}}}()
  for i ∈ 1:n
    push!(nU, 2^dim)
    U0 = tavg_pca[i,:]
    push!(U0, pcp_pca[i,1])
    push!(U, create_box(U0, collect(1:dim),tmin_pca[i,:],tmax_pca[i,:]))
  end

  return Data_p_center("nmr", n, dim, nU, U, K);

end