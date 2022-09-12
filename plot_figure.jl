using GraphRecipes, Plots, JSON

res = []
push!(res, readdlm("res/solution_p619_4_0.2_2_worst.txt"))

for data in res
    g = SimpleGraph(n)
    tree = []
    nodes = []
    
end

for i in 1:n
    set_props!(mgr, i, Dict(
        Symbol("xNode") => xNode[i],
        Symbol("yNode") => yNode[i]
        ))
end
graphplot(g,
    show = true,
    edgelabel=edgelabel,
    edge_label_box = true,
    names=nodelabel,
    fontsize = 2,
    curves = true,
    linecolor = :darkgrey,
    #nodeshape=:circle,
    nodesize=0,
    linewidth=0.2,
    x=xNode,
    y=yNode)
savefig("example_$(length(unique(values(d)))).pdf")
gplot(g, nodelabel=nodelabel)