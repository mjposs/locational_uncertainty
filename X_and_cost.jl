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
   @constraint(model, sum(x[e] for e in E) ≤ data.n-1) #avoid having too many edges in the first iterations. Does not prevent cyckes tough. To be improved.
   return model
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

   return maximum(value), u
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
   return value, u
end
