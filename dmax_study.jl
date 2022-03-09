using JuMP
using GLPK
using Gurobi

"""
    enum_solutions(U,i,y)

    Recursive function that enumerates all the possible solutions of the robust problem

 - `U` : uncertainty set 
 - `i` : current vertex under consideration
 - `y` : solution being recursively built
"""
function enum_solutions(U,i,y)
    if (i == length(U))
        res = Array{Int,2}(undef,length(U),0);
        for k in U[i]
            res = [res [y;k]];
        end
        return res;
    else
        res = Array{Int,2}(undef,length(U),0);
        for k in U[i]
            res = [res enum_solutions(U,i+1,[y;k])];
        end
        return res;
    end
end

"""
    min_ratio_exact(n,F)

    Search for the uncertainty set maximizing the approximation factor for a given edge list. 
    The uncertainty set of each vertex contains as many points as its degree, because the worst ratio can always be obtained for such instance.    With this choice, we can arbitrarily set which points of the uncertainty sets achieve maximum distance between the extremities of each eadge. This gives rise to a continuous model.   

 - `n`   : number of vertices of the graph
 - `F` : edge list
"""
function min_ratio_exact(n,F)
    degree
    # initialize the adjacency list
    A = [Int[] for i in 1:n]
    ranks = Dict();
    for (i,j) in F
        push!(A[i],j);
        push!(A[j],i);
        ranks[(i,j)] = (length(A[i]),length(A[j]));
    end
    @info "Adjacency list:  $A"

    # initialize the uncertainty sets with one point per adjacent vertex and an extra point for the optimum
    U = [Int[] for i in 1:n];
    nU = 0; # total number of points in U
    for i in 1:n
        for j in A[i]
            nU += 1;
            push!(U[i], nU);
        end
        #nU += 1;
        #push!(U[i],nU);
    end
    @info "Uncertainty sets: $U"

    indU = 1:nU; # indices of the points in U
    E = [(k,l) for (i,j) in F for k in U[i] for l in U[j]];

    ##
    model = Model(Gurobi.Optimizer);

    #
    @variable(model, dintra[i in 1:n] >= 0); 
    @variable(model, d[k in indU, l in indU] >= 0);
    @variable(model, dmax[i in 1:n, j in 1:n ; (i,j) in F] >= 0);
    @variable(model, ω >= 0);

    #
    @objective(model, Min, ω);

    # we restrict the search to uncertainty sets where the points are equidistant
    for i in 1:n
        for k in U[i], l in U[i]
            if k < l 
                @constraint(model, d[k,l] == dintra[i])
            end
        end
    end

    # Distance constraints
    # identity
    @constraint(model, ct_identity[k in indU], d[k,k] == 0);
    # symmetry 
    @constraint(model, ct_symmetry[k in indU, l in indU ; k <l], d[k,l] == d[l,k]);
    # transitivity
    for i in indU, j in indU, k in indU
        if i < j & j < k
            @constraint(model, d[i,j] + d[j,k] >= d[i,k]);
            @constraint(model, d[i,j] + d[i,k] >= d[j,k]);
            @constraint(model, d[i,k] + d[j,k] >= d[i,j]);
        end
    end

    # compute the pairwise maximum distances
    @constraint(model, ctsumdmax, sum(dmax[i,j] for (i,j) in F) == length(F));
    for (i,j) in F
        if i < j
            @constraint(model, dmax[i,j] == d[U[i][ranks[(i,j)][1]],U[j][ranks[(i,j)][2]]]);
            for k in U[i], l in U[j]
                @constraint(model, dmax[i,j] >= d[k,l]);
            end
        end
    end

    #Enumerate all the possible solutions of the robust problem
    Y = enum_solutions(U, 1, Array{Int,1}());
    for y in eachcol(Y)
        edges = [];
        for (i,j) in F
            edges = [edges ; (y[i],y[j])];
        end

        @constraint(model, ω >= sum(d[k,l] for (k,l) in edges));
    end
    optimize!(model);


    return objective_value(model)/(length(F)), value.(d), value.(dmax);
end

"""
    min_ratio_dual(n,F)

    Search for the uncertainty set maximizing the approximation factor for a given edge list. 
    Here, we consider a relaxation of the robust problem to dualize it and avoid the enumeration of an exponential number of constraints in the epigraphic formulation of the objective. 
    The function returns an upper bound on the worst-case ratio

 - `n`   : number of vertices of the graph
 - `F` : edge list
"""
function min_ratio_dual(n,F)

    # initialize the adjacency list
    A = [Int[] for i in 1:n]
    ranks = Dict();
    for (i,j) in F
        push!(A[i],j);
        push!(A[j],i);
        ranks[(i,j)] = (length(A[i]),length(A[j]));
    end
    @debug "Adjacency list:  $A"

    # initialize the uncertainty sets with one point per adjacent vertex and an extra point for the optimum
    U = [Int[] for i in 1:n];
    nU = 0; # total number of points in U
	ii = [] # table i[k]
    for i in 1:n
        for j in A[i]
            nU += 1;
            push!(U[i], nU);
			push!(ii, i)
		end
        # nU += 1;
        # push!(U[i],nU);
        # push!(ii, i)
    end
    @debug "Uncertainty sets: $U"

    indU = 1:nU; # indices of the points in U
    E = [(k,l) for (i,j) in F for k in U[i] for l in U[j]];

    ##
    model = Model(Gurobi.Optimizer);
    #
    #
    @variable(model, dintra[i in 1:n] >= 0); 
    @variable(model, d[k in indU, l in indU] >= 0);
    @variable(model, dmax[i in 1:n, j in 1:n ; (i,j) in F] >= 0);
    @variable(model, ω >= 0);

    #
    @objective(model, Min, ω);

    # we restrict the search to uncertainty sets where the points are equidistant
    for i in 1:n
        for k in U[i], l in U[i]
            if k < l 
                @constraint(model, d[k,l] == dintra[i])
            end
        end
    end

    # Distance constraints
    # identity
    @constraint(model, ct_identity[k in indU], d[k,k] == 0);
    # symmetry 
    @constraint(model, ct_symmetry[k in indU, l in indU ; k <l], d[k,l] == d[l,k]);
    # transitivity
    for i in indU, j in indU, k in indU
        if i < j & j < k
            @constraint(model, d[i,j] + d[j,k] >= d[i,k]);
            @constraint(model, d[i,j] + d[i,k] >= d[j,k]);
            @constraint(model, d[i,k] + d[j,k] >= d[i,j]);
        end
    end

    # compute the pairwise maximum distances
    @constraint(model, ctsumdmax, sum(dmax[i,j] for (i,j) in F) == length(F));
    for (i,j) in F
        if i < j
            @constraint(model, dmax[i,j] == d[U[i][ranks[(i,j)][1]],U[j][ranks[(i,j)][2]]]);
            for k in U[i], l in U[j]
                @constraint(model, dmax[i,j] >= d[k,l]);
            end
        end
    end

	#add the dual constraints
	@variable(model, α[i in 1:n])
	@variable(model, β[k in indU, j in 1:n] ≥ 0)
	@variable(model, γ[k in indU, j in 1:n] ≥ 0)
	@constraint(model, [i in 1:n, k in U[i]], α[i] - sum(β[k,j] for (i,j) in F if i==ii[k]) - sum(γ[k,j] for (j,i) in F if i==ii[k]) ≥ 0)
	@constraint(model, [(k,l) in E], β[k,ii[l]] + γ[l,ii[k]] ≥ d[k,l])
	@objective(model, Min, sum(α[i] for i in 1:n))

	#
    optimize!(model);
	@info "Dual optimal value is $(objective_value(model)/length(F))"

	# verification
    @debug begin
	   "---- running verification"
	   check = Model(Gurobi.Optimizer)
	   @variable(check, δ[k in indU], Bin)
	   @variable(check, Δ[(k,l) in E], Bin)
	   @constraint(check, [i in 1:n], sum(δ[k] for k in U[i]) == 1)
	   @constraint(check, [(i,j) in F, k in U[i]], sum(Δ[(k,l)] for (k′,l) in E if k′ == k && in(l,U[j])) ≤ δ[k])
	   @constraint(check, [(j,i) in F, k in U[i]], sum(Δ[(l,k)] for (l,k′) in E if k′ == k && in(l,U[j])) ≤ δ[k])
	   @objective(check, Max, sum(value(d[k,l])*Δ[(k,l)] for (k,l) in E))
	   optimize!(check);
	   "Check optimal value is $(objective_value(check)/length(F))"
    end

    return objective_value(model)/(length(F)), value.(d), value.(dmax);
end


######################################
# 3-PATH
######################################

# V  = [1 ; 2 ; 3];
# F = [(1,2) ; (2,3)];
# U = [1 2 ; 3 4 ; 5 6]
# rho_3path = maxratio(F,U);

# println("maximum factor for the 3-path: ", 1.0/rho_3path);


######################################
# TRIANGLE
######################################

# V  = [1 ; 2 ; 3];
# F = [(1,2) ; (1,3) ; (2,3)];
# U = [1 2 3 4 ; 5 6 7 8 ; 9 10 11 12]
# rho_triangle = maxratio(F,U);

# println("maximum factor for the triangle: ", 1.0/rho_triangle);


######################################
# 4-CYCLE
######################################

# V  = [1 ; 2 ; 3 ; 4];
# F = [(1,2) ; (2,3) ; (3,4) ; (1,4)];
# U = [1 2 ; 3 4 ; 5 6 ; 7 8]
# rho_4cycle = maxratio(F,U);

# println("maximum factor for the 4-cycle: ", 1.0/rho_4cycle);

######################################
# 5-CYCLE
######################################

# F = [(1,2) ; (2,3) ; (3,4) ; (4,5) ; (1,5)];
# U = [1 2 ; 3 4 ; 5 6 ; 7 8 ; 9 10]
# rho_5cycle = maxratio(F,U);

# println("maximum factor for the 5-cycle: ", 1.0/rho_4cycle);

# ######################################
# # 4-CLIQUE
# ######################################

# V  = [1 ; 2 ; 3 ; 4];
# F = [(1,2) ; (1,3) ; (1,4) ; (2,3) ; (2,4) ; (3,4)];
# U = [1 2 ; 3 4 ; 5 6 ; 7 8]
# rho_4clique = maxratio(F,U);

# println("maximum factor for the 4-clique: ", 1.0/rho_4cycle);


# ######################################
# # 4-STAR
# ######################################

# V  = [1 ; 2 ; 3 ; 4];
# F = [(1,2) ; (1,3) ; (1,4)];
# U = [1 2 ; 3 4 ; 5 6 ; 7 8]
# rho_4star = maxratio(F,U);

# println("maximum factor for the 4-star: ", 1.0/rho_4star);

######################################
# 7-TREE
######################################
#F = [(1,2) ; (1,3) ; (1,4) ; (2,5) ; (3,6) ; (4,7)];
#U = [1 2 3; 4 5 6 ; 7 8 9 ; 10 11 12 ; 13 14 15 ; 16 17 18 ; 19 20 21];
#
#rho_7bintree = maxratio(F,U);
#
#println("maximum factor for the 7-binary tree: ", 1.0/rho_7bintree);

######################################
# 7-BINARY-TREE
######################################
#F = [(1,2) ; (1,3) ; (2,4) ; (2,5) ; (3,6) ; (3,7)];
#U = [1 2 3; 4 5 6 ; 7 8 9 ; 10 11 12 ; 13 14 15 ; 16 17 18 ; 19 20 21];
#
#rho_7bintree = maxratio(F,U);
#
#println("maximum factor for the 7-binary tree: ", 1.0/rho_7bintree);

######################################
# 11-TREE
######################################
#F = [(1,2) ; (1,3) ; (2,4) ; (2,5) ; (3,6) ; (3,7) ; (4,8) ; (5,9) ; (6,10) ; (7,11)];
#U = [1 2 ; 3 4 ; 5 6 ; 7 8 ; 9 10 ; 11 12 ; 13 14 ; 15 16 ; 17 18 ; 19 20 ; 21 22];
#
#rho_11tree = maxratio(F,U);
#println("maximum factor for the 7-binary tree: ", 1.0/rho_11tree);


"""
test_trees()

Search the worst-case ratio for a set of (small) tree graphs
"""
function test_trees()

    @info "FOUR LEVEL BINARY TREE"
    F = [(1,2) ; (1,3) ; (2,4) ; (2,5) ; (3,6) ; (3,7) ; (4,8) ; (4,9) ; (5,10) ; (5,11) ; (6,12) ; (6,13) ; (7,14) ; (7,15)];

    @info "Exact method"
    rho_4lv_bintree = min_ratio_exact(15,F);
    @info "maximum factor for the 4-level binary tree: $(1.0/rho_4lv_bintree)"

    @info "Approximation with dualized model";
    rho_4lv_bintree_bound = min_ratio_dual(15,F);
    @info "dual bound maximum factor for the 4-level ternary tree: $(1.0/rho_4lv_bintree_bound)"


    @info "TREE LEVEL TERNARY TREE"
    F = Any[(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (2, 7), (3, 8), (3, 9), (3, 10), (4, 11), (4, 12), (4, 13)];

    @info "Exact method"
    rho_3lv_tertree = min_ratio_exact(13,F);
    @info "maximum factor for the 3-level ternary tree: $(1.0/rho_3lv_tertree)"

    @info "Approximation with dualized model";
    rho_3lv_tertree_bound = min_ratio_dual(13,F);
    @info "dual bound maximum factor for the 4-level ternary tree: $(1.0/rho_3lv_tertree_bound)"

    @info "FOUR LEVEL TERNARY TREE"
    for i in 5:13
        for j in 1:3
            global F = [F;(i, 13 + 3*(i-5) + j)]
        end
    end

    @info "Exact method"
    rho_4lv_tertree = min_ratio_exact(40,F);
    @info "maximum factor for the 4-level ternary tree: $(1.0/rho_4lv_tertree)"

    @info "Approximation with dualized model";
    rho_4lv_tertree_bound = min_ratio_dual(40,F);
    @info "dual bound maximum factor for the 4-level ternary tree: $(1.0/rho_4lv_tertree_bound)"
end


"""
test_clique(K)

Search the worst-case ratio of a K-clique

`K` : size of the largest clique that will be tested
"""
function test_clique(K)

    @info "Worst case ratio of a $(K)-CLIQUE"

    F = [];
    for i in 1:K
        for j in i+1:K
            F = [F ; (i,j)];
        end
    end


    @debug begin
        "Exact method"
        rho_exact, d, dmax = min_ratio_exact(K,F);
        "maximum factor for the $(K)-clique: $(1.0/rho_exact)"
        "\n"
        "internal distances"
        for i in 1:K
            "- in vertex $i"
            for k in (i-1)*(K-1)+1:i*(K-1)
                for l in (i-1)*(K-1)+1:i*(K-1)
                    if k < l
                        @info "d[$k,$l]=$(d[k,l])"
                    end
                end
            end
        end
        "\n"
        "pairwise distances"
        for (i,j) in F
            "- edge ($i,$j)"
            for k in (i-1)*(K-1)+1:i*(K-1)
                for l in (j-1)*(K-1)+1:j*(K-1)
                    "d[$k,$l]=$(d[k,l])"
                end
            end
        end
    end


    @info "Approximation with dualized model";
    rho_dual,d,dmax = min_ratio_dual(K,F);
    @info "dual bound maximum factor for the $(K)-clique: $(1.0/rho_dual)"
end

"""
test_quasiclique(K)

Search the worst-case ratio of a K-quasi-clique, i.e., a clique minus an edge

`K` : size of the largest clique that will be tested
"""
function test_quasiclique(K)

    @info "Worst case ratio of a $(K)-QUASI-CLIQUE"

    F = [(1,j) for j in 2:K-1];
    for i in 2:K
        for j in i+1:K
            F = [F ; (i,j)];
        end
    end

    @info "Approximation with dualized model";
    rho_dual,d,dmax = min_ratio_dual(K,F);
    @info "dual bound maximum factor for the $(K)-clique: $(1.0/rho_dual)"
end

######################################
# ENUMERATE ALL GRAPHS WITH N VERTICES
######################################
"""
enum_connected_graphs(n,(i,j),F, cur_deg, prev_deg)

Recursive function that enumerate the edge lists of all connected graphs with n vertices
Break symmetries by imposing non-increasing degrees of the vertices
"""
function enum_connected_graphs(n,(i,j), F)
    if (i==n-1)
        res = Set();
        G = Graph(n);
        for (k,l) in F
            add_edge!(G, k, l);
        end
        degrees = Graphs.degree(G)[end:-1:1];
        if !issorted(degrees) return res; end
        if is_connected(G)
            push!(res, F);
            if degrees[1] < degrees[2] push!(res, [F;(n-1,n)]); end
        else
            if degrees[1] >= degrees[2] return res; end
            add_edge!(G, n-1, n);
            if is_connected(G)
                push!(res, [F;(n-1,n)]);
            end
        end
        return res;
    else
        res = Set();
        if j == n
            union!(res, enum_connected_graphs(n,(i+1,i+2), F));
            union!(res, enum_connected_graphs(n, (i+1,i+2), [F;(i,j)]));
        else
            union!(res, enum_connected_graphs(n,(i,j+1), F));
            union!(res, enum_connected_graphs(n, (i,j+1), [F;(i,j)]));
        end
        return res;
    end
end


"""
get_all_ratios(n,k)

Get the ratios for all connected graphs with n vertices and at most k distinct points per uncertainty set
"""
function get_all_ratios(n)
    all_graphs = enum_connected_graphs(n,(1,2),[]);
    ratios = [];
    for F in all_graphs
        println("\n---------------------------");
        println("Get the ratio of graph:");
        println("\t F = ", F);
        rho = min_ratio_exact(n,F);
        push!(ratios, 1.0/rho);
        println("maximum factor : ", 1.0/rho);
    end
    return ratios;
end

function get_all_ratios_adjlist(n)
    all_graphs = enum_connected_graphs(n,(1,2),[]);
    ratios = [];
    for F in all_graphs
        println("\n---------------------------");
        println("Get the ratio of graph:");
        println("\t F = ", F);
        rho = min_ratio_exact(n,F);
        push!(ratios, 1.0/rho);
        println("maximum factor : ", 1.0/rho);
    end
    return ratios;
end
