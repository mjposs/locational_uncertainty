function create_model(data::Data, output_flag::Int64 = 0)
   model = Model()
   if Optimizer == "Gurobi"
      set_optimizer(model, () -> Gurobi.Optimizer(GUROBI_ENV))
      set_optimizer_attribute(model, "OutputFlag", output_flag)
      set_optimizer_attribute(model, "TuneOutput", 1)
      set_optimizer_attribute(model, "TimeLimit", TIMELIMIT)
      set_optimizer_attribute(model, "NodefileStart", 0.5)
   else
      set_optimizer(model, CPLEX.Optimizer)
      set_optimizer_attribute(model,"CPX_PARAM_SCRIND", output_flag)
      set_optimizer_attribute(model,"CPXPARAM_ScreenOutput", output_flag)
      #set_optimizer_attribute(model,"CPX_PARAM_MIPDISPLAY", 1)
      set_optimizer_attribute(model,"CPX_PARAM_TILIM", TIMELIMIT)
      set_optimizer_attribute(model,"CPX_PARAM_NODEFILEIND", 2)
      set_optimizer_attribute(model,"CPX_PARAM_THREADS", THREADS)
   end
   return model
end

#-------------------------------------------------------------------------------

function build_IP_model(data::Data_STP)
   E = data.E
   A = 1:2*data.m
   T = 1:data.t′
   V = 1:data.n
   model = create_model(data, 1);
   @variable(model, x[E], Bin)
   @variable(model, f[A,T] ≥ 0)
   @variable(model, ω ≥ 0)

   @objective(model, Min, ω)
   @constraint(model, [t in T, i in V], sum(f[a,t] for a in data.δ⁺[i]) - sum(f[a,t] for a in data.δ⁻[i]) == data.b[i,t])
   @constraint(model, [e in 1:length(E), t in T], f[e,t] + f[e+data.m,t] ≤ x[E[e]])
   @constraint(model, sum(x[e] for e in E) ≤ data.n-1) #avoid having too many edges in the first iterations. Does not prevent cycles though. To be improved.
   return model
end

#-------------------------------------------------------------------------------

function c(x_val, data::Data_STP)

   function fill_u(node,index)
      u[node] = index
      if indmax[node,index] != -1
         for (i,k) in indmax[node,index]
            fill_u(i,k)
         end
      end
   end

   g = SimpleGraph(data.n)
   for e in data.E x_val[e] > 1-ϵ && add_edge!(g, e[1], e[2]) end
   value = zeros(nv(g),data.nU) #value function
   indmax = convert(Matrix{Any},fill(-1,nv(g),data.nU)) #records bests
   "DP loop"
   root = data.T[1]
   dfs = dfs_tree(g, root)
   order = reverse(topological_sort_by_dfs(dfs))
   order = order[ max.( LightGraphs.degree(g)[order].>1 , order.==root ) ] # disregard leafs and isolated vertices (but keep the root)
   for i in order, k in 1:data.nU
      indmax[i,k] = []
      for j in outneighbors(dfs,i)
         index = argmax( [ value[j,l] + data.d[i,k,j,l] for l in 1:data.nU ] )
         value[i,k] += value[j,index] + data.d[i,k,j,index]
         push!(indmax[i,k],(j,index))
      end
   end
   "get the optimal solution"
   index = argmax(value[order[end],:])
   maxval = value[order[end],index]
   u = fill(-1,nv(g))
   fill_u(order[end],index)

   # some nodes may not be covered by current solution
   u = abs.(u) # arbitrarily select the first point in U[i] for each node that is not covered by x

   return maximum(value), u
end


#-------------------------------------------------------------------------------

function build_IP_model(data::Data_UFLP)
   I = data.I
   J = data.J
   model = create_model(data, 1);
   @variable(model, y[J], Bin)
   @variable(model, x[(i,j) in [(i,j) for i in I for j in J]], Bin)
   @variable(model, ω ≥ 0)

   @objective(model, Min, ω)
   @constraint(model, [i in I], sum(x[(i,j)] for j in J) == 1)
   @constraint(model, [i in I, j in J], x[(i,j)] ≤ y[j])
   @constraint(model, sum(y[j] for j in J) ≤ data.p)

   return model
end


#-------------------------------------------------------------------------------

function c(x_val, y_val, data::Data_UFLP)
   I = data.I
   J = data.J
   J′ = []
   d = data.d
   value = 0
   u = Vector{Int64}(undef,data.n)
   for j in J
      if y_val[j] < ϵ
         u[j] = 1 # take arbitrarily the first position for locations not taken
      else
         push!(J′,j)
      end
   end
   for j in J′
      I′ = [] #contains the nodes connected to J[j] in the solution
      for i in I
         x_val[(i,j)] > 1-ϵ && push!(I′,i)
      end
      best = 0
      best_pos_for_j = -1
      u_best = Vector{Int64}(undef,length(I′))
      u′ = Vector{Int64}(undef,length(I′))
      for u_j in 1:length(data.U[j])
         sum = 0
         # For each position u_j of the root of the star (the facility), we take the worst-positions for each leaf (client) I′[k]
         for k in 1:length(I′)
            u′[k] = argmax(d[j,u_j,I′[k],1:length(data.U[I′[k]])])
            sum += d[j,u_j,I′[k],u′[k]]
         end
         if sum > best
            best_pos_for_j = u_j
            best = sum
            u_best = copy(u′)
         end
      end
      u[j] = best_pos_for_j
      u[I′] = u_best
      value += best
   end

   # some nodes may not be covered by current solution
   u = abs.(u) # arbitrarily select the first point in U[i] for each node that is not covered by x
   return value, u
end



"""
   build_IP_model(data::Data_clustering)
"""
function build_IP_model(data::Data_clustering)
   n = data.n;
   K = data.K;
   V = 1:n;
   E = data.E
   n_per_cluster = floor(Int64, n/K);
   model = create_model(data, 1);
   # x_ij=1 iff i and j are in the same cluster
   @variable(model, x[(i,j) in E, Bin)
   # z_ii=1 iff vertex i is the leader of a cluster (meaning it has smallest index in the cluster); z_ij=1 iff j is in the cluster lead by i  
   @variable(model, z[i in V, j in V ; i <= j], Bin)
   # objective value in the epigraphic formulation
   @variable(model, ω ≥ 0)

   # we minimize the sum pairwise squared distance in each cluster
   @objective(model, Min, ω)

   # # each vertex must be assigned to exactly one leader or be the leader of one lcuster
   # @constraint(model, [j in V], sum(z[i,j] for i in 1:j) == 1);
   # # each vertex can be assigned only to a leader
   # @constraint(model, [i in V, j in V; j > i], z[i,j] ≤ z[i,i]);
   # # there must be K leaders
   # @constraint(model, ct_nclusters, sum(z[i,i] for i in V) == K);
   # # each cluster contains exactly K vertices
   # @constraint(model, [i in V], sum(z[i,j] for j in i+1:data.n) == (n_per_cluster - 1) * z[i,i]);
   # # i and j are in the same cluster iff they are associated to the same leader
   # @constraint(model, [h in V, i in V, j in V ; i < j && h ≤ i], x[i,j] ≥ z[h,i] + z[h,j] - 1);
   # #strengthen the formulation with valid cuts
   # # we need at least k+1 cluster to cover k*n_per_cluster+1 vertices
   # @constraint(model, [k in 1:K-1], sum(z[i,i] for i in 1:k*n_per_cluster+1) >= k+1);
   # # if j is associated to leader i, it is in the same cluster as i
   # @constraint(model, [i in V, j in V ; i < j], x[i,j] ≥ z[i,j]);
   # each vertex is with at least n_per_cluster-1 other vertices
   @constraint(model, [i in V], sum(x[(i,j)] for j in i+1:n) + sum(x[(j,i)] for j in 1:i-1) == n_per_cluster - 1);
   # transitivity constraints on x
   @constraint(model, [(i,k) in E, j in V ; j > i && j < k], x[(i,k)] ≥ x[(i,j)] + x[(j,k)] - 1);
   @constraint(model, [(i,k) in E, j in V ; j > k], x[(i,k)] ≥ x[(i,j)] + x[(k,j)] - 1);
   @constraint(model, [(i,k) in E, j in V ; j < i], x[(i,k)] ≥ x[(j,i)] + x[(j,k)] - 1);

   return model
end



#-------------------------------------------------------------------------------

function c(g::SimpleGraph{Int64}, data::Data_clustering)
   n = data.n;
   nU = data.nU;
   d² = data.d²;
   u = Int64.(zeros(n));
   total_obj_value = 0.0;
   # the solution graph provides a partition of the vertices into cliques whose worst-case costs can be computed independently
   for V in connected_components(g)
      if (length(V) != floor(Int64, n/data.K))
         @error "The connected component is not the proper size: $(length(V))"
      end
      separation = create_model(data);
      # variables indicating which position is chosen for each vertex among its uncertainty set
      @variable(separation, y[i in V, k in 1:nU[i]], Bin);
      # linearization variables: =1 for i, j, ki, kj iff position ki is chosen for vertex i and position kj is chosen for wertex j
      @variable(separation, Y[i in V, j in V, ki in 1:nU[i], kj in 1:nU[j]; i < j], Bin);

      # one position must be chosen for each vertex
      @constraint(separation, [i in V], sum(y[i,k] for k in 1:data.nU[i]) == 1);
      # pairwise indicators must be greater than each single indicator
      @constraint(separation, [i in V, j in V, ki in 1:nU[i], kj in 1:nU[j] ; i < j], Y[i,j,ki,kj] ≤ y[i,ki]);
      @constraint(separation, [i in V, j in V, ki in 1:nU[i], kj in 1:nU[j] ; i < j], Y[i,j,ki,kj] ≤ y[j,kj]);
      @constraint(separation, ct_is_edge[i in V, j in V ; i < j], sum(Y[i,j,ki,kj] for ki in 1:nU[i], kj in 1:nU[j]) ≤ 1);

      # we maximize the sum of distances between chosen positions
      @objective(separation, Max, sum(Y[i,j,ki,kj] * d²[i,j][ki,kj] for i in V, j in V, ki in 1:nU[i], kj in 1:nU[j] if i < j));

      # solve the separation problem and store the solution
      optimize!(separation);
      y_val = value.(y);
      for i in V, k in 1:nU[i]
         if y_val[i,k] >= 1 - ϵ
            u[i] = k;
         end
      end
      total_obj_value += objective_value(separation)
   end
   if (!isempty(findall(u .== 0)))
      @error "the cluster graph does not cover every vertex: $(u)"
      u[findall(u .== 0)] .= 1;
   end

   return total_obj_value, u
end

"""
   build_separation
"""
function build_separation(data::Data_clustering)
   n = data.n;
   V = 1:n
   nU = data.nU;
   d² = data.d²;
   u = Int64.(zeros(n));
   total_obj_value = 0.0;
   # the solution graph provides a partition of the vertices into cliques whose worst-case costs can be computed independently
   separation = create_model(data);
   # variables indicating which position is chosen for each vertex among its uncertainty set
   @variable(separation, y[i in V, k in 1:nU[i]], Bin);
   # linearization variables: =1 for i, j, ki, kj iff position ki is chosen for vertex i and position kj is chosen for wertex j
   @variable(separation, Y[i in V, j in V, ki in 1:nU[i], kj in 1:nU[j]; i < j], Bin);

   # one position must be chosen for each vertex
   @constraint(separation, [i in V], sum(y[i,k] for k in 1:data.nU[i]) == 1);
   # pairwise indicators must be greater than each single indicator
   @constraint(separation, [i in V, j in V, ki in 1:nU[i], kj in 1:nU[j] ; i < j], Y[i,j,ki,kj] ≤ y[i,ki]);
   @constraint(separation, [i in V, j in V, ki in 1:nU[i], kj in 1:nU[j] ; i < j], Y[i,j,ki,kj] ≤ y[j,kj]);
   @constraint(separation, ct_is_edge[i in V, j in V ; i < j], sum(Y[i,j,ki,kj] for ki in 1:nU[i], kj in 1:nU[j]) ≤ 1);

   # we maximize the sum of distances between chosen positions
   @objective(separation, Max, sum(Y[i,j,ki,kj] * d²[i,j][ki,kj] for i in V, j in V, ki in 1:nU[i], kj in 1:nU[j] if i < j));

   return separation;
end

function c(g::SimpleGraph{Int64}, data::Data_interval)
   n = data.n;
   nU = data.nU;
   u_val = zeros(n);
   total_obj_value = 0.0;
   D = 1:data.dimension;
   set_optimizer_attribute(m, "NonConvex", 2);
   # the solution graph provides a partition of the vertices into cliques whose worst-case costs can be computed independently
   for V in connected_components(g)
      if (length(V) != floor(Int64, n/data.K))
         @error "The connected component is not the proper size: $(length(V))"
      end
      separation = create_model(data);
      # variables indicating which position is chosen for each vertex among its uncertainty set
      @variable(separation, data.lb[i,d] <= u[i in V, d in D] <= data.ub[i,d]);

      @objective(separation, Max, sum((u[i,d] - u[j,d])^2 for i in V, j in V, d in D if i < j))

      # solve the separation problem and store the solution
      optimize!(separation);
      for i in V, d in D
         u_val[i,d] = value(u[i,d]);
      end
      total_obj_value += objective_value(separation)
   end
   if (!isempty(findall(u_val .== 0)))
      @error "the cluster graph does not cover every vertex: $(u)"
      u_val[findall(u_val .== 0)] .= lb[findall(u_val .== 0)];
   end

   return total_obj_value, u_val
end
