#------

function printres(time, res, seed, algo, data::Data_STP, folder::String)
    file = "res/Steiner/$(folder)$(algo).txt"
    f = open(file,"a")
    #if isfile(file)   
    #else
    #    f = open(file,"w")
    #end
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.Δ) $(data.nU) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end

function printres(time, res, seed, algo, data::Data_SPL, folder::String)
    f = open("res/SPL/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.n) $(data.m) $(length(data.I)) $(data.nU) $(length(data.I)) $(data.p) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end

function printres(time, res, seed, algo, data::Data_clustering, folder::String)
    f = open("res/clustering/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.K) $(data.n) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end


function printres(time, res, seed, algo, data::Data_p_center, folder::String)
    f = open("res/pcenter/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.K) $(data.n) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end


function compare_all_methods(data::Data,seed,folder::String)
    time_exact = @elapsed res_exact = exact(data)
    seed > 0 && printres(time_exact, res_exact, seed, "exact", data, folder)
    for costs in [("worst", data.c_max), ("center",data.c_center), ("avg",data.c_avg)]
        @info "Solve deterministic counterpart using $(costs[1])"
        time = @elapsed res = heuristic_deterministic(data, costs[2])
        seed > 0 && printres(time, res, seed, costs[1], data, folder)
    end
    if typeof(data) == Data_STP
        time = @elapsed res = heuristic_adr(data)
        seed > 0 && printres(time, res, seed, "cons", data, folder)
        # Only execute the compact formulation on the small instances
        if folder != "P6E/"
            time = @elapsed res = solve_STP_compact(data)
            seed > 0 && printres(time, res, seed, "compact", data, folder)
        end
    end
end

#------

function run_steiner()
    #folder = "small/"
	folder = "P6E/"
    #INSTANCES = readdir("data/Steiner/"*folder)
    INSTANCES =["p619.stp","p620.stp","p621.stp"]

    @warn "warming up ..."
    data = read_data_STP("data/Steiner/small.stp", 0, 1)
    compare_all_methods(data,0,folder)
    @warn "now looping over all problems"

	"NOTE: Instances small correspond to the instances format_i from the paper"
    
    for instance in INSTANCES, nU in [4,8,12], Δ in [0.2,0.4,0.6], seed in 1:5
    #for nU in [4,8,12], Δ in [0.2,0.4,0.6], seed in 1:10, size in 1:3
        Random.seed!(seed)
        #data = create_small_STP(size,Δ,nU)
        data = read_data_STP("data/Steiner/"*folder*instance, Δ, nU)
        compare_all_methods(data,seed,folder)

        #flush(stdout)
    end
end

#------

function run_SPL()
	@info "Warming up ..."
	n, m, cardI, nU, p, seed = 10, 15, 1, 1, 1, 0
    Random.seed!(0)
	data = build_SPL(n, m, cardI, nU, p)
    compare_all_methods(data,seed,"")
	@info "... finished!"

	n, cardI, p = 80, 15, 3
	maxU = Int64(floor(n/cardI))
	for n in [100,150,200], cardI in [round(n/6),round(n/5),round(n/4)], m in [n, round(n*1.5), round(n*2)], nU in [3,4,5], p in [3,4,5], seed in 1:2
        Random.seed!(seed)
		@warn "($n, $m, $cardI, $nU, $seed, $p)"
		data = build_SPL(n, m, cardI, nU, p)

        compare_all_methods(data,seed,"")
	end
end

function run_synthetic_clustering(n_per_cluster::Int64)
    U = read_cars()
    data_cars =  Data_clustering("cars", length(U), 2, [length(U[i]) for i in 1:length(U)], U, 4);
    compare_all_methods(data_cars)

    # Synthetic datasets of De Carvalho and Lechevallier (2009)
    # They propose two configurations
    # configuration 1:
    μ1_1 = [28 ; 62 ; 50 ; 57];
    μ2_1  = [23 ; 30; 15 ; 48];
    var1_1 = [144 ; 81 ; 49 ; 16];
    var2_1 = [16 ; 49 ; 81 ; 144];
    cov_1 = [0.8 ; 0.7 ; 0.6 ; 0.9];
    # configuration 2:
    μ1_2 = [28 ; 62 ; 50 ; 57];
    μ2_2  = [23 ; 30; 15 ; 37];
    var1_2 = [100 ; 81 ; 100 ; 81];
    var2_2 = [9 ; 16 ; 16 ; 9];
    cov_2 = [0.7 ; 0.8 ; 0.7 ; 0.8];

    # in the publication, they test 100 times with each of interval_length among 1, 4, 9, 14 and 19; we did not observe much variability, so we test only once for each interval length
    for width in [1.0 4.0 9.0 14.0 19.0]
        # configuration 1:
        data_config1 = build_gaussian_clustering("synth_cfg1_width=$width", n_per_cluster, width, μ1_1, μ2_1, var1_1, var2_1, cov_1);

        # configuration 2:
        data_config2 = build_gaussian_clustering("synth_cfg2_width=$width", n_per_cluster, width, μ1_2, μ2_2, var1_2, var2_2, cov_2);

        compare_all_methods(data_config1)
        compare_all_methods(data_config2)
    end

    # configuration 3:
    # μ1 = [50 ; 50];
    # μ2  = [30 ; 30];
    # data_synthetic_config3 = build_together_clustering(n_per_cluster, uncertainty_width, μ1, μ2);
end

function run_interval()
    U = read_cars()
    data_cars =  Data_p_center("cars", length(U), 2, [length(U[i]) for i in 1:length(U)], U, 4);
    compare_all_methods(data_cars)

    data_nmr = read_meteo_data()
    compare_all_methods(data_nmr)
end

function run_all()
    run_steiner();
    run_SPL();
    run_synthetic_clustering(6)
    run_interval();
end

#------

run_SPL()
# run_steiner()
