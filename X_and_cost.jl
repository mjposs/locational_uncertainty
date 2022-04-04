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
      set_optimizer_attribute(model, "CPX_PARAM_SCRIND", output_flag)
      set_optimizer_attribute(model, "CPXPARAM_ScreenOutput", output_flag)
      #set_optimizer_attribute(model,"CPX_PARAM_MIPDISPLAY", 1)
      set_optimizer_attribute(model, "CPX_PARAM_TILIM", TIMELIMIT)
      set_optimizer_attribute(model, "CPX_PARAM_NODEFILEIND", 2)
      set_optimizer_attribute(model, "CPX_PARAM_THREADS", THREADS)
   end
   return model
end

#-------------------------------------------------------------------------------

function build_IP_model(data::Data_STP)
   n = data.n
   E = data.E
   A = 1:2*data.m
   T = 1:data.t′
   V = 1:data.n
   model = create_model(data, 0)
   @variable(model, x[E], Bin)  # ∀e∈E, true if edge e is taken in the tree 
   @variable(model, f[A, T] ≥ 0)
   @variable(model, ω ≥ 0)

   @objective(model, Min, ω)
   @constraint(model, [t in T, i in V], sum(f[a, t] for a in data.δ⁺[i]) - sum(f[a, t] for a in data.δ⁻[i]) == data.b[i, t])
   @constraint(model, [e in 1:length(E), t in T], f[e, t] + f[e+data.m, t] ≤ x[E[e]])
   @constraint(model, sum(x[e] for e in E) ≤ n - 1)  # avoid having too many edges in the first iterations. Does not prevent cycles though. To be improved.
   @constraint(model, ω ≥ sum(data.c_avg[e] * x[e] for e ∈ E))
   return model
end

#-------------------------------------------------------------------------------

function build_IP_model_compact(data::Data_STP)
   n = data.n
   E = data.E
   A = 1:2*data.m
   T = 1:data.t′
   V = 1:data.n
   model = create_model(data, 0)
   @variable(model, x[E], Bin)  # ∀e∈E, true if edge e is taken in the tree 
   @variable(model, f[A, T] ≥ 0)
   @variable(model, z[V, V, 1:data.nU] ≥ 0) # z[i,k,level] worst-case cost of DP at label (i,k) at level 

   @objective(model, Min, ω)
   @constraint(model, [level in V, k in 1:data.nU], ω ≥ z[data.t′,k,level])
   @constraint(model, [i in V, level in setdiff(V,n), k in 1:data.nU], ω[i,k,level] ≥ sum(x[(i,j)]) )
   @constraint(model, [t in T, i in V], sum(f[a, t] for a in data.δ⁺[i]) - sum(f[a, t] for a in data.δ⁻[i]) == data.b[i, t])
   @constraint(model, [e in 1:length(E), t in T], f[e, t] + f[e+data.m, t] ≤ x[E[e]])
   @constraint(model, sum(x[e] for e in E) ≤ n - 1) #avoid having too many edges in the first iterations. Does not prevent cycles though. To be improved.
   @constraint(model, ω ≥ sum(data.c_avg[e] * x[e] for e ∈ E))
   return model
end

#-------------------------------------------------------------------------------

function c(x_val, data::Data_STP)
   function fill_u(node, index)
      u[node] = index
      if isassigned(indmax, node, index)
         for (i, k) in indmax[node, index]
            fill_u(i, k)
         end
      end
   end

   # create the graph induced by current solution
   n = data.n
   E = data.E
   Ex = E[findall([x_val[e] for e in E] .> 1 - ϵ)]
   nU = data.nU
   g = SimpleGraph(n)
   for e in Ex
      add_edge!(g, e[1], e[2])
   end
   Vx = findall(Graphs.degree(g) .>= 1)

   # Compute the true cost of the solution by dynamic programming
   value = zeros(nv(g), data.nU) #value function
   indmax = Matrix{Vector{Tuple{Int,Int}}}(undef, nv(g), data.nU) #records bests
   "DP loop"
   ## recover the rooted directed tree corresponding to the solution
   root = data.T[1]
   dfs = dfs_tree(g, root)
   order = reverse(topological_sort_by_dfs(dfs))
   order = order[max.(Graphs.degree(g)[order] .> 1, order .== root)] # disregard leafs and isolated vertices (but keep the root)
   for i in order
      for k in 1:data.nU
         indmax[i, k] = Vector{Tupleà)={Int,Int}}()
         for j in outneighbors(dfs, i)
            val_max, l_max = findmax(value[j, :] .+ data.cost[i, j][k, :])
            value[i, k] += val_max
            push!(indmax[i, k], (j, l_max))
         end
      end
   end
   "get the optimal solution"
   val_max, k_max = findmax(value[root, :])
   u = fill(-1, nv(g))
   fill_u(root, k_max)

   # some nodes may not be covered by current solution
   u = abs.(u) # arbitrarily select the first point in U[i] for each node that is not covered by x
   @debug begin
      println("x_val= $Ex")
      println("DP value = $(val_max)")
      println("DP value verif = $(sum(data.cost[i,j][u[i],u[j]] for (i,j) in Ex))")
      println("u_DP = $(u[Vx])")
   end


   @debug begin
      # verify the objective value by solving the separation IP
      m = Model(CPLEX.Optimizer)
      @variable(m, y[Vx, 1:nU], Bin)
      @variable(m, z[(i, j) ∈ Ex, k ∈ 1:nU, l ∈ 1:nU], Bin)
      @constraint(m, [i in Vx], sum(y[i, :]) == 1)
      @constraint(m, [(i, j) ∈ Ex], sum(z[(i, j), :, :]) == 1)
      @constraint(m, [(i, j) ∈ Ex, k ∈ 1:nU, l ∈ 1:nU], z[(i, j), k, l] ≤ y[i, k])
      @constraint(m, [(i, j) ∈ Ex, k ∈ 1:nU, l ∈ 1:nU], z[(i, j), k, l] ≤ y[j, l])
      @objective(m, Max, sum(data.cost[i, j][k, l] * z[(i, j), k, l] for (i, j) ∈ Ex, k ∈ 1:nU, l ∈ 1:nU))
      set_optimizer_attribute(m, "CPX_PARAM_SCRIND", 0)
      optimize!(m)
      println("ip value = $(objective_value(m))")
   end

   return maximum(value), u
end


#-------------------------------------------------------------------------------

function build_IP_model(data::Data_UFLP)
   I = data.I
   J = data.J
   nU = data.nU

   model = create_model(data, 1)
   @variable(model, y[J], Bin)
   @variable(model, x[(i, j) in [(i, j) for i in I for j in J]], Bin)
   @variable(model, z[J] ≥ 0)
   @variable(model, ω ≥ 0)

   @objective(model, Min, ω)
   @constraint(model, [i in I], sum(x[(i, j)] for j in J) == 1)
   @constraint(model, [i in I, j in J], x[(i, j)] ≤ y[j])
   @constraint(model, sum(y[j] for j in J) ≤ data.p)
   @constraint(model, ω ≥ sum(z))
   @constraint(model, [j in J, k in 1:nU], z[j] ≥ sum(data.c_max_from_J[j,i][k] * x[(i,j)] for i in I))

   return model
end


#-------------------------------------------------------------------------------

function c(x_val, y_val, data::Data_UFLP)
   I = data.I
   J = data.J
   Jx = []
   cost = data.cost
   value = 0
   u = Vector{Int64}(undef, data.n)
   for j in J
      if y_val[j] < ϵ
         u[j] = 1 # take arbitrarily the first position for locations not taken
      else
         push!(Jx, j)
      end
   end
   for j in Jx
      # get the nodes connected to j in the solution
      Ix = Vector{Int}()  # contains the nodes connected to J[j] in the solution
      for i in I
         x_val[(i, j)] > 1 - ϵ && push!(Ix, i)
      end

      best = 0
      best_pos_for_j = -1
      u_best = Vector{Int64}(undef, length(Ix))
      u_max = Vector{Int64}(undef, length(Ix))
      for u_j in 1:length(data.U[j])
         sum = 0
         # For each position u_j of the root of the star (the facility), we take the worst-positions for each leaf (client) Ix[k]
         for ind in 1:length(Ix)
            cost_max,u_max[ind] = findmax(cost[j, Ix[ind]][u_j, :])
            sum += cost_max
         end
         if sum > best
            best_pos_for_j = u_j
            best = sum
            u_best = copy(u_max)
         end
      end
      u[j] = best_pos_for_j
      u[Ix] = u_best
      value += best
   end

   # some nodes may not be covered by current solution
   u = abs.(u) # arbitrarily select the first point in U[i] for each node that is not covered by x
   return value, u
end


#-------------------------------------------------------------------------------

function build_IP_model(data::Data_p_center)
   E = data.E
   n = data.n
   V = 1:n
   nU = data.nU

   model = create_model(data, 1)
   @variable(model, x[(i, j) in E], Bin)
   @variable(model, z[V] ≥ 0)
   @variable(model, ω ≥ 0)

   @objective(model, Min, ω)
   @constraint(model, sum(x[(i,i)] for i in V) ≤ data.K)
   @constraint(model, [j in V], sum(x[(i, j)] for i in V) == 1)
   @constraint(model, [(i,j) in E], x[(i, j)] ≤ x[(i,i)])
   @constraint(model, ω ≥ sum(z))
   @constraint(model, [i in V, k in 1:nU[i]], z[i] ≥ sum(data.c_max_from[i,j][k] * x[(i,j)] for j in V))

   return model
end


#-------------------------------------------------------------------------------

function c(x_val, data::Data_p_center)
   n = data.n
   V = 1:n
   centers = []
   cost = data.cost
   value = 0
   u = Vector{Int64}(undef, data.n)
   for i in V
      if x_val[(i,i)] > 1 - ϵ
         push!(centers, i)
      end
   end
   @debug println("centers = $centers")
   for i in centers
      # get the nodes connected to j in the solution
      Jx = Vector{Int}()  # contains the nodes connected to centers[i] in the solution
      for j in V
         j != i && x_val[(i, j)] > 1 - ϵ && push!(Jx, j)
      end
      @debug println("Nodes in center $i : $Jx")

      best = 0
      best_pos_for_center = -1
      u_best = Vector{Int64}(undef, length(Jx))
      u_max = Vector{Int64}(undef, length(Jx))
      for u_i in 1:length(data.U[i])
         sum = 0
         # For each position u_j of the root of the star (the facility), we take the worst-positions for each leaf (client) Ix[k]
         for ind in 1:length(Jx)
            cost_max,u_max[ind] = findmax(cost[i, Jx[ind]][u_i, :])
            sum += cost_max
         end
         if sum > best
            best_pos_for_center = u_i
            best = sum
            u_best = copy(u_max)
         end
      end
      u[i] = best_pos_for_center
      u[Jx] = u_best
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
   n = data.n
   K = data.K
   V = 1:n
   E = data.E
   n_per_cluster = floor(Int64, n / K)
   model = create_model(data, 1)
   # x_ij=1 iff i and j are in the same cluster
   @variable(model, x[(i, j) in E], Bin)
   # objective value in the epigraphic formulation
   @variable(model, ω ≥ 0)

   # we minimize the sum pairwise squared distance in each cluster
   @objective(model, Min, ω)

   # each vertex is with at least n_per_cluster-1 other vertices
   @constraint(model, [i in V], sum(x[(i, j)] for j in i+1:n) + sum(x[(j, i)] for j in 1:i-1) ≥ n_per_cluster - 1)
   @constraint(model, [i in V], sum(x[(i, j)] for j in i+1:n) + sum(x[(j, i)] for j in 1:i-1) ≤ n_per_cluster)

   # transitivity constraints on x
   @constraint(model, [(i, k) in E, j in V; j > i && j < k], x[(i, k)] ≥ x[(i, j)] + x[(j, k)] - 1)
   @constraint(model, [(i, k) in E, j in V; j > k], x[(i, k)] ≥ x[(i, j)] + x[(k, j)] - 1)
   @constraint(model, [(i, k) in E, j in V; j < i], x[(i, k)] ≥ x[(j, i)] + x[(j, k)] - 1)
   @constraint(model, ω ≥ sum(data.c_avg[i,j] * x[(i,j)] for (i,j) ∈ E))

   return model
end



#-------------------------------------------------------------------------------

function c(g::SimpleGraph{Int64}, data::Data_clustering)
   n = data.n
   nU = data.nU
   cost = data.cost
   u = Int64.(zeros(n))
   total_obj_value = 0.0
   # the solution graph provides a partition of the vertices into cliques whose worst-case costs can be computed independently
   for V in connected_components(g)
      @debug println("V = $V")
      if (length(V) < floor(Int64, n / data.K))
         @error "The connected component is not the proper size: $(length(V))"
      end
      separation = build_separation(data, V)

      # solve the separation problem and store the solution
      optimize!(separation)
      y_val = value.(separation[:y])
      Y_val = value.(separation[:Y])
      for i in V, k in 1:nU[i]
         if y_val[i, k] >= 1 - ϵ
            u[i] = k
         end
      end
      total_obj_value += objective_value(separation)
   end
   if (!isempty(findall(u .== 0)))
      @error "the cluster graph does not cover every vertex: $(u)"
      u[findall(u .== 0)] .= 1
   end

   return total_obj_value, u
end

"""
   build_separation
"""
function build_separation(data::Data_clustering, V::Vector{Int})
   n = data.n
   nU = data.nU
   cost = data.cost
   u = Int64.(zeros(n))
   total_obj_value = 0.0
   # the solution graph provides a partition of the vertices into cliques whose worst-case costs can be computed independently
   separation = create_model(data, 0)
   # variables indicating which position is chosen for each vertex among its uncertainty set
   @variable(separation, y[i in V, k in 1:nU[i]], Bin)
   # linearization variables: =1 for i, j, ki, kj iff position ki is chosen for vertex i and position kj is chosen for wertex j
   @variable(separation, Y[i in V, j in V, ki in 1:nU[i], kj in 1:nU[j]; i < j], Bin)

   # one position must be chosen for each vertex
   @constraint(separation, [i in V], sum(y[i, k] for k in 1:nU[i]) == 1)
   # pairwise indicators must be greater than each single indicator
   @constraint(separation, [i in V, j in V, ki in 1:nU[i] ; i < j], sum(Y[i, j, ki, kj] for kj in 1:nU[j]) ≤ y[i, ki])
   @constraint(separation, [i in V, j in V, kj in 1:nU[j] ; i < j], sum(Y[i, j, ki, kj] for ki in 1:nU[i]) ≤ y[j, kj])
   @constraint(separation, ct_is_edge[i in V, j in V; i < j], sum(Y[i, j, ki, kj] for ki in 1:nU[i], kj in 1:nU[j]) == 1)

   # we maximize the sum of distances between chosen positions
   @objective(separation, Max, sum(Y[i, j, ki, kj] * cost[i, j][ki, kj] for i in V, j in V, ki in 1:nU[i], kj in 1:nU[j] if i < j))

   return separation
end

function c(g::SimpleGraph{Int64}, data::Data_interval)
   n = data.n
   u_val = zeros(n)
   total_obj_value = 0.0
   D = 1:data.dimension
   set_optimizer_attribute(m, "NonConvex", 2)
   # the solution graph provides a partition of the vertices into cliques whose worst-case costs can be computed independently
   for V in connected_components(g)
      if (length(V) != floor(Int64, n / data.K))
         @error "The connected component is not the proper size: $(length(V))"
      end
      separation = create_model(data)
      # variables indicating which position is chosen for each vertex among its uncertainty set
      @variable(separation, data.lb[i, dim] <= u[i in V, dim in D] <= data.ub[i, dim])

      @objective(separation, Max, sum((u[i, dim] - u[j, dim])^2 for i in V, j in V, dim in D if i < j))

      # solve the separation problem and store the solution
      optimize!(separation)
      for i in V, dim in D
         u_val[i, dim] = value(u[i, dim])
      end
      total_obj_value += objective_value(separation)
   end
   if (!isempty(findall(u_val .== 0)))
      @error "the cluster graph does not cover every vertex: $(u)"
      u_val[findall(u_val .== 0)] .= lb[findall(u_val .== 0)]
   end

   return total_obj_value, u_val
end
