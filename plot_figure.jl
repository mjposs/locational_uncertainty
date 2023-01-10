using DelimitedFiles, Graphs, Cairo, Compose, Fontconfig, GraphPlot, Statistics, Colors

res = []
nU = 4
common_nodes = [3, 28, 38, 48, 76, 66, 56, 46, 13]
basic_color = RGBA(0.7,0.7,0.7,0.5) #RGBA(0.0,0.8,0.8,1)

for algo in ["exact","center"]
    data = readdlm("res/solution_p620_4_0.6_5_$algo.txt")
    first_line_tree = findfirst(isequal(0.0), data[:,3])
    nodes_indexes = 1:(first_line_tree - 1)
    edge_indexes = first_line_tree:size(data)[1]
    nodes = Int.(data[nodes_indexes,1])
    pos = data[nodes_indexes,3:4]
    radii = data[nodes_indexes,2]
    tree = Int.(data[edge_indexes,1:2])
    degrees = [length(findall(tree .== i)) for i in nodes]
    internal_nodes = findall(degrees .> 1)
    @info mean(2*radii[internal_nodes])
    n = length(nodes)
    N = (nU+1)*n
    g = SimpleGraph(N)
    xNode = Vector{Float64}()
    yNode = Vector{Float64}()
    nodelabel = []
    nodefillc = []
    borders = []
    colors = distinguishable_colors(n)

    reverse_nodes = Dict()
    for i in 1:n
        reverse_nodes[nodes[i]] = i+(i-1)*nU
    end
    edges_color = Dict()
    edges_width = Dict()
    for e in 1:size(tree)[1]
        add_edge!(g,reverse_nodes[tree[e,1]],reverse_nodes[tree[e,2]])
        edges_color[Set((reverse_nodes[tree[e,1]],reverse_nodes[tree[e,2]]))] = RGBA(0.7,0.7,0.7,1)
        edges_width[Set((reverse_nodes[tree[e,1]],reverse_nodes[tree[e,2]]))] = 1.5
    end

    for i in 1:n
        push!(xNode,pos[i,1])
        push!(yNode,pos[i,2])
        push!(nodelabel,nodes[i])
        if in(nodes[i],common_nodes)
            push!(nodefillc, basic_color)
        else
            push!(nodefillc,colors[i])
        end
        radius = radii[i]
        if degrees[i] == 1 || in(nodes[i],common_nodes)
            radius = 0
        end
        for k in 1:nU
            push!(xNode,pos[i,1]+radius*cos(2π*k/nU))
            push!(yNode,pos[i,2]+radius*sin(2π*k/nU))
            push!(nodelabel,"")
            if in(nodes[i],common_nodes)
                push!(nodefillc, basic_color)
            else
                push!(nodefillc,colors[i])
            end
            if !in(nodes[i],common_nodes)
                if k == 1
                    origin = length(nodelabel) + nU -1
                else
                    origin = length(nodelabel) - 1
                end
                destination = length(nodelabel)
                add_edge!(g, origin, destination)
                edges_color[Set((origin,destination))] = colors[i]
                edges_width[Set((origin,destination))] = 0.25
            end
        end
    end
    
    # This seems needed as the edges are stored in an order that seems random to me
    edgestrokec = []
    edgelinewidth = []
    edgelist = collect(edges(g))
    for e in edgelist
        push!(edgestrokec, edges_color[Set((src(e),dst(e)))])
        push!(edgelinewidth, edges_width[Set((src(e),dst(e)))])
    end

    draw(PDF("$algo.pdf", 16cm, 16cm), gplot(g, xNode, yNode, nodelabel=nodelabel, nodefillc=nodefillc, nodesize=0.5, edgestrokec=edgestrokec, edgelinewidth=edgelinewidth, EDGELINEWIDTH=1.5))    
end