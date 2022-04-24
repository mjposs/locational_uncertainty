"""
exact(data)

Cutting-plane algorithm to solve the problem exactly.
"""

function exact(data)
   # build the initial relaxed IP where there is no cut to compute the worst-case cost of the solution tree
   model = build_IP_model(data)
   n = data.n
   E = data.E

   # callback function that finds the worst-case cost of the current solution and builds the corresponding "lazy" cut
   function callback(cb_data)
      status = callback_node_status(cb_data, model)
      if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
         # `callback_value(cb_data, x)` is not integer (to some tolerance).
         # If, for example, your lazy constraint generator requires an
         # integer-feasible primal solution, you can add a `return` here.
         return
      end
      # get the actual cost of the current solution and the associated worst-case position in each uncertainty set
      sep_val, u = 0, []
      x_val = Dict()
      for e in E
         x_val[e] = callback_value(cb_data, model[:x][e])
      end
      ω_val = callback_value(cb_data, model[:ω])
      if typeof(data) == Data_STP 
         "STP"

         # get the worst-case position of all the vertices by solving the separation problem
         sep_val, u = c(x_val, data)
      elseif typeof(data) == Data_SPL
         "UFLP"
         y_val = Dict()
         for j in data.J
            y_val[j] = callback_value(cb_data, model[:y][j])
         end

         # get the worst-case position of all the vertices by solving the separation problem
         sep_val, u = c(x_val, y_val, data)
      elseif typeof(data) == Data_p_center
         # get the worst-case position of all the vertices by solving the separation problem
         sep_val, u = c(x_val, data)
      elseif typeof(data) == Data_clustering
         # build the graph encoded by the current solution of the problem
         g = SimpleGraph(n)
         for e in E
            if (x_val[e] >= 1 - ϵ)
               add_edge!(g, e[1], e[2])
            end
         end

         # get the worst-case position of all the vertices by solving the separation problem
         sep_val, u = c(g, data)
      end

      @debug begin
         println("Separation value = $(sep_val), ω = $(ω_val)")
         Ex = E[findall([x_val[e] for e in E] .> 1-ϵ)]
         println("x_val = $(Ex)")
         for e in Ex
            println("cost[$(e[1]), $(e[2])] = $(data.cost[e[1], e[2]][u[e[1]], u[e[2]]] )")
         end
         sep_val_verif = sum(data.cost[i, j][u[i], u[j]] * x_val[(i, j)] for (i,j) in Ex)
         println("sep_val_verif = $sep_val_verif")
      end
      # add the optimality cut if it is not satisfied by the current solution
      if sep_val > ω_val + ϵ
         # build the corresponding optimality cut
         con = @build_constraint(model[:ω] ≥ sum(data.cost[i, j][u[i], u[j]] * model[:x][(i, j)] for (i, j) in E))
         MOI.submit(model, MOI.LazyConstraint(cb_data), con)
         ncuts += 1
         # print something every 100 cuts generated
         rem(ncuts, 100) == 0 && "$ncuts cuts generated"
      end
   end

   # Don't use lazy cut callback for the compact formulation of the SPL model
   typeof(data) != Data_SPL && MOI.set(model, MOI.LazyConstraintCallback(), callback)
   ncuts = 0

   optimize!(model)
   
   @info "Solution found of cost $(value(model[:ω])), generating $ncuts cuts"
   for e in E
      value.(model[:x][e]) .> ϵ && @debug "Edge from $(e[1]) to $(e[2])"
   end
   return objective_value(model), ncuts, MOI.get(model, MOI.RelativeGap())
end

#-------------------------------------------------------------------------------

"""
heuristic_deterministic

Classical approach ignoring uncertainty and taking a given deterministic cost for each edge. This cost will typically reflect etiher the maximum, the average or the distance between the centers of the uncertainty sets.
"""

function heuristic_deterministic(data, cost::Dict{Tuple{Int64,Int64},Float64})
   typeof(data) == Data_STP && @info "deterministic solution for $(data.instance) with Δ=$(data.Δ) and $(data.nU) extreme points"
   typeof(data) == Data_SPL && @info "deterministic solution for $(data.instance) with |I|=$(length(data.I)), |J|=$(length(data.J)), p=$(data.p)"
   typeof(data) == Data_clustering && @info "deterministic solution for $(data.instance) with  $(data.nU) extreme points"

   # initialize the model
   model = build_IP_model(data)
   E = data.E
   n = data.n

   # only the objective is modified in the heuristic solving a deterministic model
   @objective(model, Min, sum(cost[e] * model[:x][e] for e in E))
   optimize!(model)
   @info "Solution found with heuristic cost $(objective_value(model))"

   # get the true cost of the solution
   x_val = Dict()
   for e in E
      x_val[e] = value(model[:x][e])
   end
   if typeof(data) == Data_STP
      "STP"
      truecost, u = c(x_val, data)
   elseif typeof(data) == Data_SPL
      "UFLP"
      y_val = Dict()
      for j in data.J
         y_val[j] = value(model[:y][j])
      end
      truecost, u = c(x_val, y_val, data)
   elseif typeof(data) == Data_p_center
      for i in 1:n
         for j in 1:i
            x_val[(i,j)] = value(model[:x][(i,j)])
         end
      end
      truecost, u = c(x_val, data)
   elseif typeof(data) == Data_clustering
      # compute the true cost of the solution
      g = SimpleGraph(n)
      for e in E
         if x_val[e] > ϵ
            add_edge!(g, e[1], e[2])
         end
      end
      truecost, u = c(g, data)
   else
      @error("Unknow data type")
   end

   @info "True cost is $truecost"
   for e in E
      x_val[e] > ϵ && @debug "$(e[1]) - $(e[2])"
   end
   return truecost, objective_value(model), MOI.get(model, MOI.RelativeGap())
end

#-------------------------------------------------------------------------------

"""
ADR based approach from de Ruiter et al. Only for Steiner Tree
"""

function heuristic_adr(data::Data_STP)
   @info "adr heuristic for $(data.instance) with Δ=$(data.Δ) and $(data.nU) extreme points"
   V = 1:data.n
   E = data.E
   # In the following two definitions, we take only the edges with indexes less than m, meaning a single directed edge
   # instead of two opposite ones.
   δ⁺ = [ E[data.δ⁺[i][data.δ⁺[i] .≤ data.m]] for i in V ] 
   δ⁻ = [ E[data.δ⁻[i][data.δ⁻[i] .≤ data.m]] for i in V ]
   s = data.nU
   model = build_IP_model(data)
   add_bridge(model, MOI.Bridges.Constraint.SOCtoNonConvexQuadBridge)

   @variable(model, μ0[V])
   @variable(model, μ[E, 1:2]) #μ_ij
   @variable(model, νfrom[E, 1:s] ≥ 0) # used in the epigraphic reformulation of ||x_ij u_i^k - μ_ij||_2
   @variable(model, νto[E, 1:s] ≥ 0) # used in the epigraphic reformulation of ||x_ji u_i^k - μ_ji||_2

   # ν_i,ij^k ≥ ||x_ij u_i^k - μ_ij||_2
   for i in V, e in δ⁺[i], k in 1:s
      vector_for_SOCP = Vector{GenericAffExpr{Float64,VariableRef}}(undef, 3)
      vector_for_SOCP[1] = νfrom[e, k]
      vector_for_SOCP[2:3] = [model[:x][e]*data.U[i][k][q] - μ[e, q] for q in 1:2]
      @constraint(model, vector_for_SOCP in SecondOrderCone())
   end
   # ν_i,ji^k ≥ ||x_ji u_i^k - μ_ji||_2
   for i in V, e in δ⁻[i], k in 1:s
      vector_for_SOCP = Vector{GenericAffExpr{Float64,VariableRef}}(undef, 3)
      vector_for_SOCP[1] = νto[e, k]
      vector_for_SOCP[2:3] = [model[:x][e]*data.U[i][k][q] - μ[e, q] for q in 1:2]
      @constraint(model, vector_for_SOCP in SecondOrderCone())
   end
   @constraint(model, model[:ω] ≥ sum(μ0[i] for i in V))
   @constraint(model, [i in V, k in 1:s], μ0[i] ≥ sum(νfrom[e,k] for e in δ⁺[i]) + sum(νto[e,k] for e in δ⁻[i]))

   optimize!(model)
   @info "Solution found of cost $(objective_value(model))"
   x_val = Dict()
   for e in E
      x_val[e] = value(model[:x][e])
   end
   truecost, u = c(x_val, data)
   @info "True cost is $truecost"

   @debug begin
      println("x: $(value.(model[:x]))")
      println("μ0: $(value.(μ0))")
      println("μ: $(value.(μ))")
      println("νfrom: $(value.(νfrom))")
      println("νto: $(value.(νto))")
      for i in V, k in 1:s
         println("node $i, position $k: $(sum([value(νfrom[e,k]) for e in δ⁺[i]]) + sum([value(νto[e,k]) for e in δ⁻[i]]))")
      end
   end
   
   return truecost, objective_value(model), MOI.get(model, MOI.RelativeGap())
end

#-------------------------------------------------------------------------------
"""
Compact exact formulation for the Steiner Tree Problem. Relies on a directed formulation instead
of the non-directed formulation used in the build_IP_model(). Hence, rather to add undirected variables
to make the formulation compatible with the function exact() above, we consider this formulation
as a full *solution agorithm*.
"""

function solve_STP_compact(data::Data_STP)
   n = data.n
   m = data.m
   from = [data.from; data.to] # append the reverse edges
   to = [data.to; data.from] # append the reverse edges
   A = [(from[a], to[a]) for a in 1:2m]
   δ⁺ = [A[data.δ⁺[i]] for i in 1:n]
   δ⁻ = [A[data.δ⁻[i]] for i in 1:n]
   T = data.T
   T0 = T[1:data.t′]
   V = 1:data.n
   r = T[data.t]
   cardU = 1:data.nU

   model = create_model()
   @variable(model, x[A], Bin)  # ∀a∈A, true if directed edge a is taken in the arborescence
   @variable(model, f[A, T0] ≥ 0)
   @variable(model, z[V, cardU] ≥ 0) # z[i,k] worst-case cost of DP at label (i,k)
   @variable(model, Z[A, cardU] ≥ 0) # used to linearize the maximization over ℓ
   @variable(model, X[A, cardU] ≥ 0) # linearize the product x[a]*Z[A,k,ℓ]
   @variable(model, ω ≥ 0)

   @objective(model, Min, ω)
   @constraint(model, [k in cardU], ω ≥ z[r,k])
   @constraint(model, [i in V,k in cardU,ℓ in cardU], z[i,k] ≥ sum(X[a,k] for a in δ⁺[i]))
   @constraint(model, [a in A, k in cardU, ℓ in cardU], x[a] => {X[a,k] ≥ Z[a,k]}) #Notice big-M could be used instead
   @constraint(model, [a in A, k in cardU, ℓ in cardU], Z[a,k] ≥ data.cost[a[1],a[2]][k,ℓ] + z[a[2],ℓ])

   @constraint(model, [t in T0, i in V], sum(f[a, t] for a in δ⁺[i]) - sum(f[a, t] for a in δ⁻[i]) == data.b[t][i])
   @constraint(model, [a in A, t in T0], f[a, t] ≤ x[a])
   optimize!(model)

   @info "Solution found of cost $(objective_value(model))"
   
   @debug begin 
      println("Print solution:")
      println("Edges")
      for a in A
         value(x[a]) > 0.001 && println("$a")
      end
      println("Nodes values")
      for i in V, k in cardU
         value(z[i,k])>0.001 && println("node $i, position $k: $(value(z[i,k]))")
      end
   end

   return objective_value(model), -1,  MOI.get(model, MOI.RelativeGap())
end