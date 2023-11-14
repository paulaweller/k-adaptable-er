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
    termination_plot_from_csv(filename::String; terminated="both")

    filename without .csv ending, terminated is one from both, bb, box, neither, bb_infeas
"""
function termination_plot_from_csv(filename::String; terminated="both", k=2)

    # extract data from file, remove lines with strings where bb is infeas
    alldata_dirty_allk = DataFrame(CSV.File("$(filename).csv"))
    alldata_dirty = alldata_dirty_allk[alldata_dirty_allk[!,:k] .== k,:]
    alldata = alldata_dirty[alldata_dirty[!,:it_bb] .!= "s",:]
    alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
    alldata[!,:θ_bb] = parse.(Float64, alldata[!,:θ_bb])
    # instances where bb/box terminated/ didn't terminate
    bb_term, bb_unterm = split_time(alldata, "bb")
    box_term, box_unterm = split_time(alldata, "box")

    # instances where bb was infeasible
    bb_infeas = alldata_dirty[alldata_dirty[!,:it_bb] .== "s",:]

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
        obj_plot = plot(ylabel="objective", title="BB infeasible")
        θ_box_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :θ_box]
        θ_box_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :θ_box]
        θ_pb_term = plot_data[plot_data[!,:runtime_box] .<= 3600, :θ_pb]
        θ_pb_unterm = plot_data[plot_data[!,:runtime_box] .> 3600, :θ_pb]
        if length(θ_box_term) > 0
            boxplot!(obj_plot, ["Box t. $(length(θ_box_term))"], [θ_box_term], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            boxplot!(obj_plot, ["PB"], [θ_pb_term], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
        end
        if length(θ_box_unterm) > 0
            boxplot!(obj_plot, ["Box unt. $(length(θ_box_unterm))"], [θ_box_unterm], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            boxplot!(obj_plot, ["PB"], [θ_pb_unterm], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
        end
        savefig(obj_plot, "source/plots/termination/k$(k)_bb_infeasible_obj.pdf")

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
        savefig(time_plot, "source/plots/termination/k$(k)_bb_infeasible_time.pdf")
    else
        obj_plot_bb_box = plot(ylabel="objective", title="$(nrow(plot_data)) instances")
        θ_bb = plot_data[:, :θ_bb]
        θ_box = plot_data[:, :θ_box]
        θ_pb = plot_data[:, :θ_pb]
        if nrow(plot_data) > 0
            boxplot!(obj_plot_bb_box, ["BB" "Box" "PB"], [θ_bb θ_box θ_pb], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
            savefig(obj_plot_bb_box, "source/plots/termination/k$(k)_$(terminated)_terminate_obj.pdf")
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
            savefig(time_plot, "source/plots/termination/k$(k)_$(terminated)_terminate_time.pdf")
        end
    end
    # combi = plot(obj_plot_bb_box, time_plot, layout=(1,2), legend=false)
    # savefig(combi, "$(filename)_boxplot.pdf")

end

"""
    k_plot_from_csv(filename::String; method="box")

    filename without .csv ending
"""
function k_plot_from_csv(filename::String; method="box", relative=false)

    # extract data from file, remove lines with strings where bb is infeas
    alldata_dirty = DataFrame(CSV.File("$(filename).csv"))
    if method == "bb"
        # filter feasible ones
        alldata = alldata_dirty[alldata_dirty[!,:it_bb] .!= "s",:]
        alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
        alldata[!,:θ_bb] = parse.(Float64, alldata[!,:θ_bb])
        term, unterm = split_time(alldata, "bb")
        θ_key = :θ_bb
        time_key = :runtime_bb
    elseif method == "box"
        alldata = alldata_dirty
        term, unterm = split_time(alldata, "box")
        θ_key = :θ_box
        time_key = :runtime_box
    elseif method == "pb"
        term = alldata_dirty
        θ_key = :θ_pb
        time_key = :runtime_pb
    else
        error("method is invalid, must be one of bb, box, pb")
    end
    # sort by k
    term_k1 = term[term[!,:k] .== 1,:]
    term_k2 = term[term[!,:k] .== 2,:]
    term_k3 = term[term[!,:k] .== 3,:]

    term_is1 = semijoin(term_k1, term_k2, on = [:instance, :n, :m, :pc])
    term_is1 = semijoin(term_k1, term_k3, on = [:instance, :n, :m, :pc])
    sort!(term_is1, [:n, :m, :pc, :instance])

    term_is2 = semijoin(term_k2, term_k1, on = [:instance, :n, :m, :pc])
    term_is2 = semijoin(term_k2, term_k3, on = [:instance, :n, :m, :pc])
    sort!(term_is2, [:n, :m, :pc, :instance])

    term_is3 = semijoin(term_k3, term_k1, on = [:instance, :n, :m, :pc])
    term_is3 = semijoin(term_k3, term_k2, on = [:instance, :n, :m, :pc])
    sort!(term_is3, [:n, :m, :pc, :instance])

    if relative == true 
        obj_plot = plot(ylabel="relative improvement of objective")
        θ1 = 100 .*(1 .- term_is1[:, θ_key]./term_is1[:, θ_key])
        θ2 = 100 .*(1 .- term_is2[:, θ_key]./term_is1[:, θ_key])
        θ3 = 100 .*(1 .- term_is3[:, θ_key]./term_is1[:, θ_key])
    elseif relative == false
        obj_plot = plot(ylabel="objective")
        θ1 = term_is1[:, θ_key]
        θ2 = term_is2[:, θ_key]
        θ3 = term_is3[:, θ_key]
    else 
        error("value for relative is invalid")
    end


    boxplot!(obj_plot, ["k=1" "k=2" "k=3"], [θ1 θ2 θ3], leg=false, linewidth=2,linecolour= :match,fillalpha = 0.4)
    if relative == true 
        savefig(obj_plot, "source/plots/k_comparison/k3_$(method)_relobj.pdf")
    else 
        savefig(obj_plot, "source/plots/k_comparison/k3_$(method)_obj.pdf")
    end

    t1 = term_is1[:, time_key]
    t2 = term_is2[:, time_key]
    t3 = term_is3[:, time_key]
    time_plot = plot(ylabel="runtime")
    boxplot!(time_plot, ["k=1" "k=2" "k=3"], [t1 t2 t3], 
        leg=false, 
        linewidth=2,
        linecolour= :match,
        fillalpha = 0.4)

    savefig(time_plot, "source/plots/k_comparison/k3_$(method)_time.pdf")
    # combi = plot(obj_plot_bb_box, time_plot, layout=(1,2), legend=false)
    # savefig(combi, "$(filename)_boxplot.pdf")

end