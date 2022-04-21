using DelimitedFiles

names = ["\\exact","\\worst","\\hcenter","\\avg","\\cons","\\compact"]
NU = 1
MAXVAL = 60
folder = "small"
#folder = ""


#------

function print_costs_STP()
    algos = ["worst","center","avg","cons"]
    #if folder != "small_i"
    #    global algos = algos[1:2] #remove heur_cons for large instances
    #end
    datafile = readdlm("res/Steiner/$folder/exact.txt")
    D = sort(unique(datafile[:,2]))
    nseed = maximum(datafile[:,8])
    ninstances = length(unique(datafile[:,1]))
    nU = length(unique(datafile[:,3]))
    values = Vector{Vector{Float64}}(undef,length(algos))

    exact = datafile[1:end,5]
    
    for i in 1:length(algos)
        datafile = readdlm("res/Steiner/$folder/$(algos[i]).txt")
        values[i] = 100*(datafile[1:end,5] - exact)./exact
    end
    for Δ in D
        indexes = findall(datafile[:,2] .== Δ)
        val = [ values[i][indexes] for i in 1:length(algos) ]
        if length(algos) == 2
            allvalues = sort(union(val[1],val[2],MAXVAL))
        else
            allvalues = sort(union(val[1],val[2],val[3],MAXVAL))
        end
        res = Vector{Matrix{Any}}(undef,length(algos))
        for i in 1:length(algos)
            val[i] = sort(val[i])
            res[i] = ["x" "solved"]
        end
    	for value in allvalues, i in 1:length(algos)
    		number = findlast(val[i] .<= value)
    		if number !== nothing res[i] = [ res[i] ; value 100*number/(nseed*ninstances*nU) ]
    		else res[i] = [ res[i] ; value 0 ]
    		end
    	end
        for i in 1:length(algos)
            writedlm("res/Steiner/$(folder)_$(algos[i])_$(Δ).txt", res[i])
        end
    end
end


#------

function print_costs_UFLP()
    MAXVAL = 20
    algos = ["center","worst"]
    datafile = readdlm("res/UFLP/exact.txt")
    M = sort(unique(datafile[:,3]))
    nseed = maximum(datafile[:,end])
    NU = sort(unique(datafile[:,4]))
    values = Vector{Vector{Float64}}(undef,2)

    exact = datafile[1:end,7]
    for i in 1:length(algos)
        datafile = readdlm("res/UFLP/$(algos[i]).txt")
        values[i] = 100*(datafile[1:end,7] - exact)./exact
    end

    for nU in NU
        ninstances = 3
        indexes = findall(datafile[:,4] .== nU)
        val = [ values[i][indexes] for i in 1:length(algos) ]
        allvalues = sort(union(val[1],val[2],MAXVAL))
        res = Vector{Matrix{Any}}(undef,3)
        for i in 1:length(algos)
            val[i] = sort(val[i])
            res[i] = ["x" "solved"]
        end
    	for value in allvalues, i in 1:length(algos)
    		number = findlast(val[i] .<= value)
    		if number !== nothing res[i] = [ res[i] ; value 100*number/(nseed*ninstances) ]
    		else res[i] = [ res[i] ; value 0 ]
    		end
    	end
        for i in 1:length(algos)
            writedlm("res/UFLP/$(algos[i])_NU_$(nU).txt", res[i])
        end
    end

    for m in M
        ninstances = 3
        indexes = findall(datafile[:,3] .== m)
        val = [ values[i][indexes] for i in 1:length(algos) ]
        allvalues = sort(union(val[1],val[2],MAXVAL))
        res = Vector{Matrix{Any}}(undef,3)
        for i in 1:length(algos)
            val[i] = sort(val[i])
            res[i] = ["x" "solved"]
        end
        for value in allvalues, i in 1:length(algos)
            number = findlast(val[i] .<= value)
            if number !== nothing res[i] = [ res[i] ; value 100*number/(nseed*ninstances) ]
            else res[i] = [ res[i] ; value 0 ]
            end
        end
        for i in 1:length(algos)
            writedlm("res/UFLP/$(algos[i])_M_$(m).txt", res[i])
        end
    end
end

#------

function print_times_STP()
    parameters = ["Delta", "nU"]
    folder == "small" && push!(parameters, "dim")
    column = Dict()
    column["Delta"] = 2
    column["nU"] = 3
    column["dim"] = 1
    algos = ["worst","center","avg","cons","exact"]
    folder == "small" && push!(algos, "compact")
    for param in parameters
        col = column[param]
        datafile = readdlm("res/Steiner/$folder/exact.txt")
        if param == "dim"
            datafile[:,1] = parse.(Int64,[ datafile[i,1][end] for i in 1:size(datafile)[1] ])
        end
        allvalues = sort(unique(datafile[:,col]))
        table = ["x"; allvalues]
        for algo in algos
            datafile = readdlm("res/Steiner/$folder/$algo.txt")
            if param == "dim"
                datafile[:,1] = parse.(Int64,[ datafile[i,1][end] for i in 1:size(datafile)[1] ])
            end
            allvalues = sort(unique(datafile[:,col]))
            newcol = Vector{Any}(undef,1)
            newcol[1] = algo
            for val in allvalues
                rows_to_sum = findall(datafile[:,col] .== val)
                mean = sum(datafile[rows_to_sum, 4])/length(rows_to_sum)
                push!(newcol, mean)
            end
            table = [table newcol]
        end
        writedlm("res/Steiner/$(folder)_times_$(param).txt", table)
    end
end

#------

function print_times_UFLP()
    parameters = ["nU", "m"]
    column = Dict()
    column["nU"] = 4
    column["m"] = 3
    algos = ["exact","center","worst"]
    for param in parameters
        col = column[param]
        datafile = readdlm("res/UFLP/exact.txt")
        if param == parameters[1]
            datafile[:,1] = parse.(Int64,[ datafile[i,1][end] for i in 1:size(datafile)[1] ])
        end
        allvalues = sort(unique(datafile[:,col]))
        table = ["x"; allvalues]
        for algo in algos
            datafile = readdlm("res/UFLP/$algo.txt")
            if param == parameters[1]
                datafile[:,1] = parse.(Int64,[ datafile[i,1][end] for i in 1:size(datafile)[1] ])
            end
            allvalues = sort(unique(datafile[:,col]))
            newcol = Vector{Any}(undef,1)
            newcol[1] = algo
            for val in allvalues
                rows_to_sum = findall(datafile[:,col] .== val)
                mean = sum(datafile[rows_to_sum, 6])/length(rows_to_sum)
                push!(newcol, mean)
            end
            table = [table newcol]
        end
        writedlm("res/UFLP/times_$(param).txt", table)
    end
end

print_costs_STP()
print_times_STP()

#print_costs_UFLP()
#print_times_UFLP()

function print_costs()
    algos = ["heur_det","heur_dmax","heur_cons"]
    if folder != "small"
        global algos = algos[1:2] #remove heur_cons for large instances
    end
    datafile = readdlm("res/Steiner/"*folder*"/exact.txt")
    D = sort(unique(datafile[:,2]))
    nseed = maximum(datafile[:,8])
    ninstances = length(unique(datafile[:,1]))
    values = Vector{Vector{Float64}}(undef,3)

    exact = datafile[1:end,5]
    for i in 1:length(algos)
        datafile = readdlm("res/Steiner/"*folder*"/$(algos[i]).txt")
        values[i] = 100*(datafile[1:end,5] - exact)./exact
    end
    for Δ in D
        indexes = findall(datafile[:,2] .== Δ)
        val = [ values[i][indexes] for i in 1:length(algos) ]
        if length(algos) == 2
            allvalues = sort(union(val[1],val[2],MAXVAL))
        else
            allvalues = sort(union(val[1],val[2],val[3],MAXVAL))
        end
        res = Vector{Matrix{Any}}(undef,3)
        for i in 1:length(algos)
            val[i] = sort(val[i])
            res[i] = ["x" "solved"]
        end
        for value in allvalues, i in 1:length(algos)
            number = findlast(val[i] .<= value)
            if number !== nothing res[i] = [ res[i] ; value 100*number/(nseed*ninstances) ]
            else res[i] = [ res[i] ; value 0 ]
            end
        end
        for i in 1:length(algos)
            writedlm("res/Steiner/$(folder)_$(algos[i])_$(Δ).txt", res[i])
        end
    end
end

function print_times(type)
    #parameters = ["Delta", "nU", "dim"]
    parameters = ["Delta"]
    column = Dict()
    column["Delta"] = 2
    column["nU"] = 3
    column["dim"] = 1
    #algos = ["exact","heur_det","heur_dmax","heur_cons"]
    algos = ["exact","heur_det","heur_dmax"]
    for param in parameters
        col = column[param]
        datafile = readdlm("res/Steiner/$type/exact.txt")
        if param == "dim"
            datafile[:,1] = parse.(Int64,[ datafile[i,1][end] for i in 1:size(datafile)[1] ])
        end
        allvalues = sort(unique(datafile[:,col]))
        table = ["x"; allvalues]
        for algo in algos
            datafile = readdlm("res/Steiner/$type/$algo.txt")
            if param == "dim"
                datafile[:,1] = parse.(Int64,[ datafile[i,1][end] for i in 1:size(datafile)[1] ])
            end
            allvalues = sort(unique(datafile[:,col]))
            newcol = Vector{Any}(undef,1)
            newcol[1] = algo
            for val in allvalues
                rows_to_sum = findall(datafile[:,col] .== val)
                mean = sum(datafile[rows_to_sum, 4])/length(rows_to_sum)
                push!(newcol, mean)
            end
            table = [table newcol]
        end
        writedlm("res/Steiner/$(type)_times_$(param).txt", table)
    end

end

#print_costs()
#print_times("P6E")
