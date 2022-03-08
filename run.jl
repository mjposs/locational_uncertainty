#------

function printres(time, res, seed, algo, data::Data_STP)
    f = open("res/Steiner/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.Δ) $(data.nU) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end

function printres(time, res, seed, algo, data::Data_UFLP)
    f = open("res/UFLP/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.n) $(data.m) $(data.nU) $(length(data.I)) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end

function printres(time, res, seed, algo, data::Data_clustering)
    f = open("res/clustering/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.K) $(data.n) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end


function printres(time, res, seed, algo, data::Data_p_center)
    f = open("res/pcenter/$(algo).txt","a")
    if algo[1] == 'P'
        algo = algo[4:end]
    end
    s = string("$(data.instance) $(data.K) $(data.n) $(@sprintf("%.2f",time))  $(@sprintf("%.2f",res[1]))  $(@sprintf("%.2f",res[2]))  $(@sprintf("%.2f",res[3]))  $(seed)\n")
    write(f,s)
    println(s)
    close(f)
end


function compare_all_methods(data::Data)
    seed = 1;
    time_exact = @elapsed res_exact = exact(data);
    time_worst = @elapsed res_worst = heuristic_deterministic(data, data.c_max);
    time_center = @elapsed res_center = heuristic_deterministic(data, data.c_center);
    time_avg = @elapsed res_avg = heuristic_deterministic(data, data.c_avg);

    printres(time_exact, res_exact, seed, "exact", data);
    printres(time_worst, res_worst, seed, "worst", data);
    printres(time_center, res_center, seed, "center", data);
    printres(time_avg, res_avg, seed, "avg", data);

    #time_adr = @elapsed res_adr = heuristic_adr(data);
    #printres(time_adr, res_adr, seed, "adr", data);
end

#------

function run_steiner()
    Random.seed!(1)
	seed = 0

    @warn "warming up ..."
    # data = read_data_STP("data/Steiner/small.stp", 0, 1)
    # @info "CP for $(data.instance) with Δ=$(data.Δ) and $(data.nU) extreme points"
    # compare_all_methods(data)

    @warn "now looping over all problems"

	"NOTE: Instances small and small_i correspond to the instances format_i from the paper"
    #folder = "small_i/"
	folder = "P6E/"
    #INSTANCES = readdir("data/Steiner/"*folder)
    INSTANCES =["p619.stp"] #,"p620.stp","p621.stp"]
    for instance in INSTANCES, nU in [10], Δ in 0.1
    #for nU in [5,10,20], Δ in [0.1,0.5,1], seed in 1:20, size in 1:2
        #@warn "seed number $seed"
        data = read_data_STP("data/Steiner/"*folder*instance,Δ,nU)
        #data = create_small_STP(size,Δ,nU)

        compare_all_methods(data)
        flush(stdout)
    end
end

#------

function run_UFLP()
	@info "Warming up ..."
	n, m, cardI, nU, p, seed = 10, 20, 1, 1, 1, 1
	data = build_UFLP(n, m, cardI, nU, p, seed)
    compare_all_methods(data)
	@info "... finished!"

	n, cardI, p = 80, 15, 3
	maxU = Int64(floor(n/cardI))
	for seed in 1:1, m in [n, round(n*1.5), round(n*2)], nU in [3,4,5]
		@warn "($n, $m, $cardI, $nU, $seed, $p)"
		data = build_UFLP(n, m, cardI, nU, seed, p)

        compare_all_methods(data)
	end
end

function run_clustering(n_per_cluster::Int64)
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
    run_UFLP();
    run_synthetic_clustering(6);
    run_interval();
end

#------

# run_UFLP()
# run_steiner()
