using LinearAlgebra, Dates, Random, Statistics, JuMP, Plots,StatsPlots
#StatsPlots,

"""
    obs_obj(uncertainty_set, find_plan, q_all, c)

Calculate actually observable objectives.
"""
function obs_obj(uncertainty_set, find_plan, q_all, c)

    obs_objectives = zeros(length(q_all), length(uncertainty_set))

    # for every iteration
    for i in 1:length(q_all)
        count = 0

        # for every demand scenario
        for u in uncertainty_set
            count = count+1
            # find the best plan
            obs_obj, obs_q = find_plan(u, q_all[i]..., c)
            # add it to the list
            obs_objectives[i,count] = obs_obj

        end
    end
    return obs_objectives
end

"""
    k_curve(obj_values, k_number; comp_adapt=0, static=obj_values[1], observ=0, n_val=0, m_val=0, d_max=5, p=0.5, rel_val=false)

Plot the k-curve for the given objective values.
"""
function k_curve(obj_values, k_number; comp_adapt=0, static=obj_values[1], observ=0, n_val=0, m_val=0, d_max=5, p=0.5, rel_val=false)

    # for a relative value curve (objective values relative to static objective)
    if rel_val == true

        obj_values = 100*obj_values/obj_values[1]
        observ = 100*observ/obj_values[1]
        yaxes = "%"

    else

        yaxes = "objective"

    end

    # plot the static objective
    static_obj = []
    for i in 1:length(k_number)
        push!(static_obj, static)
    end
    plot(k_number, static_obj, lw = 3, color=RGB(210/255, 22/255, 53/255), label = "static", xlabel = "k", ylabel= yaxes)

    # if given, plot observable costs for each scenario in the uncertainty set
    if observ !== 0
        plot!(k_number, observ, line= (0.3, 1), color=RGB(252/255, 204/255, 95/255), label="")
        M = mean(observ, dims=2)
        plot!(k_number, M, lw=3, color=RGB(0,0,0), label="observable mean")
    end

    # plot the k-adaptable objective
    plot!(k_number, obj_values, lw = 3, color=RGB(0/255, 71/255, 119/255), label = "k-adapt")
    
    # plot the completely adaptable objective 207, 159, 205
    if comp_adapt !== 0
        comp_adapt_obj = []
        for i in 1:length(k_number)
            push!(comp_adapt_obj, comp_adapt)
        end
        plot!(k_number, comp_adapt_obj, lw = 3, color=RGB(207/255, 159/255, 205/255),label = "comp. adapt.")
    end
    # save plot to file
    savefig("results/kcurve.pdf")
end

"""
    k_curve_obs(k_number, observ; n_val=0, m_val=0, d_max=5, p=0.5, plot_scenarios=true)

Plot the k-curves of observable objectives with mean and 0.1/0.9 quantiles.
"""
function k_curve_obs(k_number, observ; n_val=0, m_val=0, d_max=5, p=0.5, plot_scenarios=true)

    # plot observable costs
    if plot_scenarios==true
        it = hcat("observable objective",fill("",1,size(observ, 2)-1))
        plot(k_number, observ,
                line= (0.5, 0.5), 
                color=RGB(122/255, 200/255, 255/255), # light blue
                #color = RGB(207/255, 159/255, 205/255), # purple
                label=it,
                xlabel = "k", ylabel= "objective")
    end
    # calculate mean of observable costs
    M = mean(observ, dims=2)

    # plot mean
    plot!(k_number, M, lw=3, color=RGB(210/255, 22/255, 53/255), label="observable mean")

    # compute and plot quantiles
    quantiles = zeros(length(k_number), 2)
    for k in 1:length(k_number)
        temp = quantile(observ[k,:], [0.1, 0.9])
        # println("temp $temp")
        quantiles[k, 1] = temp[1]
        quantiles[k, 2] = temp[2]
    end
    plot!(k_number, quantiles, lw=3, color=RGB(0/255, 71/255, 119/255), style=:dash, label=["0.1 quantile" "0.9 quantile"])
    savefig("results/kcurve_obs.pdf")
end

"""
    k_curve_multi(mult; n_val=3, m_val=6, d_max=5, p=0.5, number_of_plans=false)

Plot k-curves of multiple instances (mult = number of instances). number_of_plans indicates whether the x-axis shows the number of cells or 
the true number of plans.
"""
function k_curve_multi(mult; n_val=3, m_val=6, d_max=5, p=0.5, number_of_plans=false, rel_val=false)

    # for a relative value curve (objective values relative to static objective)
    if rel_val == true
        yaxes = "%"
    else
        yaxes = "objective"
    end

    if number_of_plans == true
        plan_plot = plot(title = "n=$n_val, m=$m_val, d_max = $d_max, p=$p", xlabel = "number of plans", ylabel= "objective")
    end
    cell_plot = plot(title = "n=$n_val, m=$m_val, d_max = $d_max, p=$p", xlabel = "number of cells", ylabel= "objective")

    # file for saving the instances in the plot
    io = open("results/kcurve_multi.txt", "w")
    write(io, "instances:\n")
    close(io)

    for i in 1:mult
        println("iteration $i")
        no_sp = n_val
        no_dp = m_val

        # randomly generate seed
        se = rand(1000:4000)

        # seed 2 was abnormally long 
        @info "seed: $se"

        # generate instance
        loc_I_gen, loc_J_gen, demand_bound, cont_perc, agg_supply_bound = generate_instance(no_sp,no_dp, se, plot_loc=false)

        io = open("results/kcurve_multi.txt", "a")
        write(io, "iteration $i: n = $no_sp, m = $no_dp, seed = $se\n")
        close(io)

        #calculate costs
        c_gen = reshape([norm(loc_I_gen[:,i]-loc_J_gen[:,j]) for j in 1:no_dp for i in 1:no_sp],no_sp,no_dp)

        # calculate k-adaptable solution
        o_v, w_v, q_v, p_v, p_true_v = k_adapt_solution(10, loc_I_gen, loc_J_gen, agg_supply_bound, demand_bound, cont_perc)

        # for a relative value curve (objective values relative to static objective)
        if rel_val == true
            o_v = 100*o_v/o_v[1]
        end

        # plot the objectives
        if number_of_plans == true 
            plot!(plan_plot, 
                p_true_v, o_v, 
                line = (0.3,3), 
                color=RGB(52/255, 89/255, 149/255), 
                label = "") 
        end 
        plot!(cell_plot,
                p_v, o_v, 
                line = (0.3,3), 
                color=RGB(52/255, 89/255, 149/255), 
                label = "") 
    end

    # save plots to file
    if number_of_plans==true
        savefig(plan_plot, "results/kcurve_multi_plans")
    end
    savefig(cell_plot, "results/kcurve_multi_cells.pdf")
end

function box_plot_from_textfiles(filenames, labelnames)

    plot_data_obj = zeros(50,length(filenames))
    plot_data_time = zeros(50, length(filenames))
    count = 1
    for o_file in filenames
        # extract objective values from file
        o_values = open(o_file) do file
            obj = []
            for ln in eachline(file)
                index = findlast(isequal('%'), ln)
                val = parse(Float64, ln[index+2:end-1])
                push!(obj, val)
            end
            obj
        end
        plot_data_obj[:, count] =  o_values
        # extract runtime values from file
        t_values = open(o_file) do file
            times = []
            for ln in eachline(file)
                index_s = findlast(isequal(':'), ln)
                index_e = findlast(isequal('s'), ln)
                val = parse(Float64, ln[index_s+2:index_e-2])
                push!(times, val)
            end
            times
        end
        plot_data_time[:, count] = t_values
        count = count+1
    end
    obj_plot = plot(xlabel="m", ylabel="objective %")
    time_plot = plot(xlabel="m", ylabel="time in s")
    boxplot!(obj_plot, labelnames, plot_data_obj, leg=false, linewidth=2,colour = [RGB(122/255, 200/255, 255/255) RGB(0/255, 71/255, 119/255) RGB(207/255, 159/255, 205/255) RGB(210/255, 22/255, 53/255)],linecolour= :match,fillalpha = 0.4)
    boxplot!(time_plot, labelnames, plot_data_time, leg=false, linewidth=2,colour = [RGB(122/255, 200/255, 255/255) RGB(0/255, 71/255, 119/255) RGB(207/255, 159/255, 205/255) RGB(210/255, 22/255, 53/255)],linecolour= :match,fillalpha = 0.4)
    savefig(obj_plot, "results/boxplot_obj.pdf")
    savefig(time_plot, "results/boxplot_time.pdf")
end

function filter_floats!(df)
    dd = filter(x -> x isa Float64, df)
    return dd
end

function split_time(df, alg)
    term = DataFrame()
    unterm = DataFrame()
    if alg=="bb"
        term = df[df[!,:runtime_bb] .<= 3600 , :]
        unterm = df[df[!,:runtime_bb] .> 3600 , :]
    end
    if alg=="box"
        term = df[df[!,:runtime_box] .<= 3600 , :]
        unterm = df[df[!,:runtime_box] .> 3600 , :]
    end
    return term, unterm
end

"""
    termination_plot_from_csv(filename::String; terminated="both", k=2, rel_to_pb=false)

    filename without .csv ending, terminated is one from both, bb, box, neither, bb_infeas
"""
function termination_plot_from_csv(filename::String; terminated="both", K=[2], rel_to_pb = false)

    # extract data from file, remove lines with strings where bb is infeas
    alldata_dirty_allk = DataFrame(CSV.File("$(filename).csv"))

    alldata_dirty = alldata_dirty_allk[in.(alldata_dirty_allk[!,:k],(K,)),:]

    alldata = alldata_dirty[alldata_dirty[!,:θ_bb] .!= 1e20,:]
    # alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
    # alldata[!,:θ_bb] = parse.(Float64, alldata[!,:θ_bb])
    # instances where bb/box terminated/ didn't terminate
    bb_term, bb_unterm = split_time(alldata, "bb")
    box_term, box_unterm = split_time(alldata, "box")

    # instances where bb was infeasible
    bb_infeas = alldata_dirty[alldata_dirty[!,:θ_bb] .== 1e20,:]

    # all combinations of termination
    both_term = semijoin(bb_term, box_term, on = [:instance, :n, :m, :pc, :k])
    only_bb_term = bb_term[bb_term[!,:runtime_box] .> 3600 , :]
    only_box_term = box_term[box_term[!,:runtime_bb] .> 3600, :]
    neither_term = bb_unterm[bb_unterm[!,:runtime_box] .> 3600, :]

        # check that no instance is missing
    if nrow(bb_infeas) + nrow(both_term) + nrow(only_bb_term) + nrow(only_box_term) + nrow(neither_term) != nrow(alldata_dirty)
        error("Missing instances: Some were not infeasible, terminated or unterminated.")
    end

    if terminated == "both"
        plot_data = both_term
    elseif terminated == "bb"
        plot_data = only_bb_term
    elseif terminated == "box"
        plot_data = only_box_term
    elseif terminated == "neither"
        plot_data = neither_term
    elseif terminated == "bb_infeas"
        plot_data = bb_infeas
    else 
        error("Option terminated is not valid. Enter terminated= one of both,bb,box,neither,bb_infeas.")
    end

    if terminated == "bb_infeas"
        
        θ_box_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :θ_box]
        θ_box_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :θ_box]
        θ_pb_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :θ_pb]
        θ_pb_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :θ_pb]

        obs_box_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :obs_box]
        obs_box_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :obs_box]
        obs_pb_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :obs_pb]
        obs_pb_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :obs_pb]
        if rel_to_pb == false
            obj_plot = plot(ylabel="objective", title="BB infeasible")
            if size(θ_box_term,1) > 0
                boxplot!(obj_plot, ["Box t. $(length(θ_box_term))" "Box t. obs"], [θ_box_term obs_box_term], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
                boxplot!(obj_plot, ["PB" "PB obs"], [θ_pb_term obs_pb_term], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            end
        elseif rel_to_pb == true
            obj_plot = plot(ylabel="relative improvement %", title="BB infeasible")
            if size(θ_box_term,1) > 0
                boxplot!(obj_plot, ["Box t. $(length(θ_box_term))" "Box t. obs"], [(100 .*(θ_pb_term.-θ_box_term)./θ_pb_term) (100 .*(obs_pb_term.-obs_box_term)./obs_pb_term)], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            end
        end
        if size(θ_box_unterm,1) > 0
            if rel_to_pb == false
                boxplot!(obj_plot, ["Box unt. $(length(θ_box_unterm))" "Box unt. obs"], [θ_box_unterm obs_box_unterm], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
                boxplot!(obj_plot, ["PB" "PB obs"], [θ_pb_unterm obs_pb_unterm], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            elseif rel_to_pb == true
                boxplot!(obj_plot, ["Box unt. $(length(θ_box_unterm))" "Box unt. obs"], [(100 .*(θ_pb_unterm.-θ_box_unterm)./θ_pb_unterm) (100 .*(obs_pb_unterm.-obs_box_unterm)./obs_pb_unterm)], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            end
        end
        rel_to_pb ? savefig(obj_plot, "source/plots/termination/k$(K)_bb_infeasible_obj_rel_to_pb.pdf") : savefig(obj_plot, "source/plots/termination/k$(K)_bb_infeasible_obj.pdf")

        time_box_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :runtime_box]
        time_box_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :runtime_box]
        time_pb_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :runtime_pb]
        time_pb_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :runtime_pb]
        time_plot = plot(ylabel="runtime", title="BB infeasible")
        if length(time_box_term) > 0
            boxplot!(time_plot, ["Box t." "PB"], [time_box_term time_pb_term], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        end
        if length(time_box_unterm) > 0
            boxplot!(time_plot, ["Box unt." "PB"], [time_box_unterm time_pb_unterm], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        end
        savefig(time_plot, "source/plots/termination/k$(K)_bb_infeasible_time.pdf")
    else
        
        θ_bb = plot_data[:, :θ_bb]
        θ_box = plot_data[:, :θ_box]
        θ_pb = plot_data[:, :θ_pb]
        obs_bb = plot_data[:, :obs_bb]
        obs_box = plot_data[:, :obs_box]
        obs_pb = plot_data[:, :obs_pb]
        if nrow(plot_data) > 0
            if rel_to_pb == false
                obj_plot_bb_box = plot(ylabel="objective", title="$(nrow(plot_data)) instances")
                boxplot!(obj_plot_bb_box, ["BB" "Box" "PB"], [θ_bb θ_box θ_pb], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4) #TODO relative
                savefig(obj_plot_bb_box, "source/plots/termination/k$(k)_$(terminated)_terminate_obj.pdf")
            elseif rel_to_pb == true
                obj_plot_bb_box = plot(ylabel="relative improvement %", title="$(nrow(plot_data)) instances")
                boxplot!(obj_plot_bb_box, ["BB" "BB obs" "Box" "Box obs"], [(100 .*(θ_pb.-θ_bb)./θ_pb) (100 .*(obs_pb.-obs_bb)./obs_pb) (100 .*(θ_pb.-θ_box)./θ_pb) (100 .*(obs_pb.-obs_box)./obs_pb)], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4) #TODO relative
                savefig(obj_plot_bb_box, "source/plots/termination/k$(K)_$(terminated)_terminate_obj_rel_obs.pdf")
            end
        end

        time_bb = plot_data[:, :runtime_bb]
        time_box = plot_data[:, :runtime_box]
        time_pb = plot_data[:, :runtime_pb]
        time_plot = plot(ylabel="runtime", title="$(nrow(plot_data)) instances")
        if nrow(plot_data) > 0
            boxplot!(time_plot, ["BB" "Box" "PB"], [time_bb time_box time_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
            savefig(time_plot, "source/plots/termination/k$(K)_$(terminated)_terminate_time.pdf")
        end
    end
    # combi = plot(obj_plot_bb_box, time_plot, layout=(1,2), legend=false)
    # savefig(combi, "$(filename)_boxplot.pdf")

end

"""
    k_plot_from_csv(filename::String; method="box")

    filename without .csv ending
"""
function k_plot_from_csv(filename::String; method="box", relative=false, observable=false)

    # extract data from file, remove lines with strings where bb is infeas
    alldata_dirty = DataFrame(CSV.File("$(filename).csv"))
    if method == "bb"
        # filter feasible ones
        alldata = alldata_dirty#[alldata_dirty[!,:it_bb] .!= "s",:]
        # alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
        # alldata[!,:θ_bb] = parse.(Float64, alldata[!,:θ_bb])
        if observable 
            term = alldata
            θ_key = :obs_bb
        else
            term, unterm = split_time(alldata, "bb")
            θ_key = :θ_bb
        end
        
        time_key = :runtime_bb
    elseif method == "box"
        alldata = alldata_dirty
        if observable
            term = alldata
            θ_key = :obs_box
        else
            term, unterm = split_time(alldata, "box")
            θ_key = :θ_box
        end
        time_key = :runtime_box
    elseif method == "pb"
        term = alldata_dirty
        observable ? (θ_key = :obs_pb) : (θ_key = :θ_pb)
        time_key = :runtime_pb
    else
        error("method is invalid, must be one of bb, box, pb")
    end
    # sort by k
    term_k1 = term[term[!,:k] .== 1,:]
    @show size(term_k1)
    term_k2 = term[term[!,:k] .== 2,:]
    @show size(term_k2)
    term_k3 = term[term[!,:k] .== 3,:]
    @show size(term_k3)

    if observable 
        term_is1 = term_k1
        term_is2 = term_k2
        term_is3 = term_k3
        sort!(term_is1, [:n, :m, :pc, :instance])
        sort!(term_is2, [:n, :m, :pc, :instance])
        sort!(term_is3, [:n, :m, :pc, :instance])
    else
        term_is1 = semijoin(term_k1, term_k2, on = [:instance, :n, :m, :pc])
        term_is1 = semijoin(term_is1, term_k3, on = [:instance, :n, :m, :pc])
        sort!(term_is1, [:n, :m, :pc, :instance])

        #term_is2 = semijoin(term_k2, term_k1, on = [:instance, :n, :m, :pc])
        term_is2 = semijoin(term_is2, term_k3, on = [:instance, :n, :m, :pc])
        sort!(term_is2, [:n, :m, :pc, :instance])

        #term_is3 = semijoin(term_k3, term_k1, on = [:instance, :n, :m, :pc])
        term_is3 = semijoin(term_is3, term_k2, on = [:instance, :n, :m, :pc])
        sort!(term_is3, [:n, :m, :pc, :instance])
    end

    if relative == true 
        obj_plot = plot(ylabel="relative improvement of objective")
        #θ1 = 100 .*(1 .- term_is1[:, θ_key]./term_is1[:, θ_key])
        θ2 = 100 .*(1 .- term_is2[:, θ_key]./term_is2[:, :θ_pb])
        θ3 = 100 .*(1 .- term_is3[:, θ_key]./term_is3[:, :θ_pb])
    elseif relative == false
        obj_plot = plot(ylabel="objective")
        #θ1 = term_is1[:, θ_key]
        θ2 = term_is2[:, θ_key]
        θ3 = term_is3[:, θ_key]
    else 
        error("value for relative is invalid")
    end
    n_inst = length(θ2)

    boxplot!(obj_plot, ["k=2" "k=3"], [θ2 θ3], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
    if relative == true 
        observable ? savefig(obj_plot, "source/plots/k_comparison/obsrelobj_k_$(method)_$(n_inst).pdf") : savefig(obj_plot, "source/plots/k_comparison/relobj_k_$(method)_$(n_inst).pdf")
    else 
        observable ? savefig(obj_plot, "source/plots/k_comparison/obsobj_k_$(method)_$(n_inst).pdf") : savefig(obj_plot, "source/plots/k_comparison/obj_k_$(method)_$(n_inst).pdf")
    end
    if observable == false 
        t1 = term_is1[:, time_key]
        t2 = term_is2[:, time_key]
        t3 = term_is3[:, time_key]
        time_plot = plot(ylabel="runtime")
        boxplot!(time_plot, ["k=1" "k=2" "k=3"], [t1 t2 t3], 
            leg=false, 
            linewidth=2,
            linecolour= :match,
            fillalpha = 0.4)

        savefig(time_plot, "source/plots/k_comparison/time_k_$(method)_$(n_inst).pdf")
    end
    # combi = plot(obj_plot_bb_box, time_plot, layout=(1,2), legend=false)
    # savefig(combi, "$(filename)_boxplot.pdf")

end

function plot_pc_vs_time(filename; time=false, objective=false)
    alldata = DataFrame(CSV.File("$(filename).csv"))

    # replace!(alldata.runtime_bb, "i" => "3600.0")
    # alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
    
    data1 = alldata[alldata[!,:pc] .== 0.1,:]
    data3 = alldata[alldata[!,:pc] .== 0.3,:]

    if time == true
        time_plot1 = plot(xlabel="pc = 0.1")
        boxplot!(time_plot1, ["bb" "box" "pb"], [data1[!,:runtime_bb] data1[!,:runtime_box] data1[!,:runtime_pb]], 
            leg=false, 
            linewidth=2,
            linecolour= :match,
            fillalpha = 0.4)

            time_plot3 = plot(xlabel="pc = 0.3")
            boxplot!(time_plot3, ["bb" "box" "pb"], [data3[!,:runtime_bb] data3[!,:runtime_box] data3[!,:runtime_pb]], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)

        cplot = plot(time_plot1, time_plot3, ylabel="runtime")
        savefig(cplot, "source/plots/percentage/time.pdf")
    end
    if objective == true
        alldata_clean = alldata[alldata[!,:θ_bb] .!= 1e20,:]
        # alldata_clean = alldata[alldata[!,:it_bb] .!= "s",:]
        # alldata_clean[!,:θ_bb] = parse.(Float64, alldata_clean[!,:θ_bb])

        data1 = alldata_clean[alldata_clean[!,:pc] .== 0.1,:]
        data3 = alldata_clean[alldata_clean[!,:pc] .== 0.3,:]

        data1 = semijoin(data1, data3, on = [:instance, :n, :m])
        data3 = semijoin(data3, data1, on = [:instance, :n, :m])

        time_plot1 = plot(xlabel="pc = 0.1")
        boxplot!(time_plot1, ["bb" "box" "pb"], [data1[!,:θ_bb] data1[!,:θ_box] data1[!,:θ_pb]], 
            leg=false, 
            linewidth=2,
            linecolour= :match,
            fillalpha = 0.4)

            time_plot3 = plot(xlabel="pc = 0.3")
            boxplot!(time_plot3, ["bb" "box" "pb"], [data3[!,:θ_bb] data3[!,:θ_box] data3[!,:θ_pb]], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)

        cplot = plot(time_plot1, time_plot3, ylabel="objective")
        savefig(cplot, "source/plots/percentage/objective.pdf")
    end
end

function plot_evol(zetas; xlimits=[0,3600], relative=true, name="zeta", last=false)
    time_plot = plot(xlims=xlimits)
    for (key, values) in zetas
        # Extract x and y coordinates from each 2D vector
        if length(values) > 1 && values != "Vector{Float64}[]"
            coordinates = vecfromstr(values)
            
            # Extract x and y coordinates from each 2D vector
            x_vals = coordinates[:,1]
            if coordinates[1,2] == 0
                if last == true
                    relative ? (y_vals = (coordinates[end,2] .- coordinates[:,2])./coordinates[end,2]) : (y_vals = coordinates[:,2])
                else 
                    relative ? (y_vals = (coordinates[2,2] .- coordinates[:,2])./coordinates[2,2]) : (y_vals = coordinates[:,2])
                end
            else 
                relative ? (y_vals = (coordinates[1,2] .- coordinates[:,2])./coordinates[1,2]) : (y_vals = coordinates[:,2])
            end
            # Plot the data
            plot!(time_plot, x_vals[1:(end)], y_vals[1:(end)], label="",line = (0.1,1), color=palette(:default)[1])
        end
    end
    if last == true
        relative ? savefig(time_plot, "source/plots/gaps/$(name)_end_rel_$(xlimits[2]).pdf") : savefig(time_plot, "source/plots/gaps/$(name)_$(xlimits[2]).pdf")
    else
        relative ? savefig(time_plot, "source/plots/gaps/$(name)_rel_$(xlimits[2]).pdf") : savefig(time_plot, "source/plots/gaps/$(name)_$(xlimits[2]).pdf")
    end
end

"""
    plot_zeta_distr_from_csv(filename::String; method="box")

    filename without .csv ending, interval= number of intervals whichin occurences are counted together?
"""
function plot_zeta_distr_from_csv(filename::String; status="terminated", n_interval=10, relative=false, K=[1,2,3])

    # extract data from file, remove lines with strings where bb is infeas
    alldata = DataFrame(CSV.File("$(filename).csv"))

    terminated, unterminated = split_time(alldata, "box")
    θ_key = :zeta_box

    if status == "terminated"
        term = terminated
    elseif status == "unterminated"
        term = unterminated
    else 
        error("status invalid, enter terminated or unterminated")
    end

    obj_plot = plot(ylabel="occurrences", xlabel="zeta")
    interval_borders = LinRange(0, 50, n_interval+1)

    # sort by k
    for k in K
        term_k = term[term[!,:k] .== k,:]

        if relative == true 
            # get first zetas, divide by them
        elseif relative == false
            θk = term_k[:, θ_key]
        else 
            error("value for relative is invalid")
        end
        
        z = []
        println(z)
        for i in 1:(n_interval)
            push!(z, number_in_interval(interval_borders[i], interval_borders[i+1], θk))
        end
        # Find nonzero entries
        lastk = findall(x -> x>0, z)

        scatter!(obj_plot, (36/(n_interval)).*lastk, z[lastk], label="k=$k")
    end
    if status =="terminated"
        relative ? savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_termin_rel_$(K).pdf") : savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_termin_$(K).pdf")
    else 
        relative ? savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_unterm_rel_$(K).pdf") : savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_untermin_$(K).pdf")
    end
end

number_in_interval(lb, rb, zetas) = count(i-> ((i <= rb) && (i > lb)), zetas)


function plot_size_vs_time(filename; percentage=0.1, K=[2])
    alldata = DataFrame(CSV.File("$(filename).csv"))

    # replace!(alldata.runtime_bb, "i" => "3600.0")
    # alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
    
    data = alldata[alldata[!,:pc] .== percentage,:]
    plots = []
    for n in [8,6,4]
        for m in [10,15,20]
            if n == 4 && m != 10
                time_plot = plot(xlabel="m = $m")#,tickfontsize=1, tickfontcolor=:white)
            elseif m ==10
                if n == 4
                    time_plot=plot(ylabel="n = $n", xlabel="m = $m")#,tickfontsize=1, tickfontcolor=:white)
                else 
                    time_plot = plot(ylabel="n = $n")#,tickfontsize=1, tickfontcolor=:white)
                end
            else
                time_plot = plot()#tickfontsize=1, tickfontcolor=:white)
            end
            # filter data
            data_n = data[data[!,:n] .== n,:]
            data_m = data_n[data_n[!,:m] .== m,:]
            data_k = data_m[in.(data_m[!,:k],(K,)),:]
            boxplot!(time_plot, ["bb" "box" "pb"], [data_k[:,:runtime_bb] data_k[:,:runtime_box] data_k[:,:runtime_pb]], 
                leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            push!(plots, time_plot)
        end
    end
    cplot = plot(plots..., layout=(3,3))
    savefig(cplot, "source/plots/size/runtime_pc$(percentage)_K$(K).pdf")
end

function observable_plot(o_pb, o_bb, o_box, oo_pb, oo_box,k, pc; rel=true)
    if rel == true
        plot1 = plot(ylabel="observable improvement")
        boxplot!(plot1, ["bb ($(length(o_bb)))" "box"], [100 .*(o_pb.-o_bb)./o_pb 100 .*(o_pb.-o_box)./o_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        if length(oo_pb) > 0
            boxplot!(plot1, ["bb infeas ($(length(oo_box)))"], [100 .*(oo_pb.-oo_box)./oo_pb], 
                    leg=false, 
                    linewidth=2,
                    linecolour= :match,
                    fillalpha = 0.4)
        end
        savefig(plot1, "source/plots/observable/rel_objective_$(pc)_k$(k).pdf")
    elseif rel == false
        plot1 = plot(ylabel="observable worst-case objective")
        boxplot!(plot1, ["bb($(length(o_bb)))" "box" "pb"], [o_bb o_box o_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        if length(oo_pb) > 0
                boxplot!(plot1, ["bb infeas ($(length(oo_box))): box" " bb infeas: pb"], [oo_box oo_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        end
        savefig(plot1, "source/plots/observable/objective_$(pc)_k$(k).pdf")
    else error("optinal argument rel must be true or false")
    end
end
