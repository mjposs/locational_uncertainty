"""
exact(data)

Cutting-plane algorithm to solve the problem exactly.
"""

function exact(data)
   # build the initial relaxed IP where there is no cut to compute the worst-case cost of the solution tree
   model = build_IP_model(data)
   n = data.n
   V = 1:n
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
      elseif typeof(data) == Data_UFLP
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

      @info "Separation value = $(sep_val), ω = $(ω_val)"
      @debug begin
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

   MOI.set(model, MOI.LazyConstraintCallback(), callback)
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
   typeof(data) == Data_UFLP && @info "deterministic solution for $(data.instance) with |I|=$(length(data.I)), |J|=$(length(data.J)), p=$(data.p)"
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
   elseif typeof(data) == Data_UFLP
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

"ADR based approach from de Ruiter et al."

function heuristic_adr(data::Data_STP)
   @info "adr heuristic for $(data.instance) with Δ=$(data.Δ) and $(data.nU) extreme points"
   V = 1:data.n
   E = data.E
   δ⁺ = [data.δ⁺[i][data.δ⁺[i].≤data.m] for i in V]
   δ⁻ = [data.δ⁻[i][data.δ⁻[i].≤data.m] for i in V]
   s = data.nU
   M = [maximum(data.cost[i, j]) for i in 1:data.n, j in 1:data.n]
   model = build_IP_model(data)
   add_bridge(model, MOI.Bridges.Constraint.SOCtoNonConvexQuadBridge)

   @variable(model, μ0[V])
   @variable(model, μ[1:data.n, 1:data.m, 1:2]) #using sparse definitions might be faster?
   @variable(model, ν[1:data.m] ≥ 0) # = ||μ_i,ij + μ_j,ij||_2
   @variable(model, νfrom[V, 1:data.m, 1:s] ≥ 0) # = ||u_i^k - μ_i,ij||_2
   @variable(model, νto[V, 1:data.m, 1:s] ≥ 0) # = ||u_i^k - μ_i,ji||_2
   @variable(model, πfrom[V, 1:data.m, 1:s] ≥ 0) # = x_e × νto[i,e,k]
   @variable(model, πto[V, 1:data.m, 1:s] ≥ 0) # = x_e × νfrom[i,e,k]

   for e in 1:data.m
      vector_for_SOCP = Vector{GenericAffExpr{Float64,VariableRef}}(undef, 3)
      vector_for_SOCP[1] = ν[e]
      vector_for_SOCP[2] = μ[data.from[e], e, 1] + μ[data.to[e], e, 1]
      vector_for_SOCP[3] = μ[data.from[e], e, 2] + μ[data.to[e], e, 2]
      @constraint(model, vector_for_SOCP in SecondOrderCone())
   end
   for i in V, e in δ⁺[i], k in 1:s
      vector_for_SOCP = Vector{GenericAffExpr{Float64,VariableRef}}(undef, 3)
      vector_for_SOCP[1] = νfrom[i, e, k]
      vector_for_SOCP[2:3] = [data.U[i][k][q] - μ[i, e, q] for q in 1:2]
      @constraint(model, vector_for_SOCP in SecondOrderCone())
   end
   for i in V, e in δ⁻[i], k in 1:s
      vector_for_SOCP = Vector{GenericAffExpr{Float64,VariableRef}}(undef, 3)
      vector_for_SOCP[1] = νto[i, e, k]
      vector_for_SOCP[2:3] = [data.U[i][k][q] + μ[i, e, q] for q in 1:2]
      @constraint(model, vector_for_SOCP in SecondOrderCone())
   end
   @constraint(model, model[:ω] ≥ sum(μ0[i] for i in V) + sum(ν[e] for e in 1:data.m))
   @constraint(model, [i in V, k in 1:s], μ0[i] ≥ sum(πfrom[i, e, k] for e in δ⁺[i]) + sum(πto[i, e, k] for e in δ⁻[i]))
   # TODO: it seems that there is an index error below, is it really E[e] and M[e] that are relevant?
   @constraint(model, [i in V, k in 1:s, e in δ⁺[i]], πfrom[i, e, k] ≥ νfrom[i, e, k] - M[e] * (1 - model[:x][E[e]]))
   @constraint(model, [i in V, k in 1:s, e in δ⁻[i]], πto[i, e, k] ≥ νto[i, e, k] - M[e] * (1 - model[:x][E[e]]))

   @debug model

   optimize!(model)
   @info "Solution found of cost $(objective_value(model))"
   x_val = Dict()
   for e in E
      x_val[e] = value(model[:x][e])
   end
   truecost, u = c(x_val, data)
   @info "True cost is $truecost"

   @debug begin
      "x $(value.(model[:x]))"
      "μ0: $(value.(μ0))"
      "μ: $(value.(μ))"
      "ν: $(value.(ν))"
      "νfrom: $(value.(νfrom))"
      "νto: $(value.(νto))"
      "πfrom: $(value.(πfrom))"
      "πto: $(value.(πto))"
   end
   return truecost, objective_value(model), MOI.get(model, MOI.RelativeGap())
end
