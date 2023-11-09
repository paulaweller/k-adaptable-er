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
    filter!(x -> x isa Float64, df)
end

function filter_time!(df)
    filter!(x-> x < 3600, df)
end

"""
    box_plot_from_csv(filename::String, labelname::String; objectives=true, times=true)

    filename without .csv ending
"""
function box_plot_from_csv(filename::String, labelname::String; objectives=true, times=true)

    plot_data_obj = zeros(50)
    plot_data_time = zeros(50)

    # extract objective values from file
    alldata = DataFrame(CSV.File("$(filename).csv"))

    if objectives == true
        θ_bb = alldata[:, :θ_bb]
        filter_floats!(θ_bb)
        if length(θ_bb) != 50
            θ_bb = zeros(50)
        end
        θ_box = alldata[:, :θ_box]
        θ_pb = alldata[:, :θ_pb]
        obj_plot = plot(xlabel="method", ylabel="objective", title=labelname)
        boxplot!(obj_plot, ["BB" "Box" "PB"], [θ_bb θ_box θ_pb], 
            leg=false, 
            linewidth=2,
            colour = [RGB(122/255, 200/255, 255/255) RGB(0/255, 71/255, 119/255) RGB(207/255, 159/255, 205/255) RGB(210/255, 22/255, 53/255)],
            linecolour= :match,
            fillalpha = 0.4)
        savefig(obj_plot, "$(filename)_boxplot.pdf")
    end
    if times==true
        runtime_bb = alldata[:, :runtime_bb]
        filter_floats!(runtime_bb)
        if length(runtime_bb) != 50
            runtime_bb = zeros(50)
        end
        runtime_box = alldata[:, :runtime_box]
        box_optimal = length(filter(x->x<3600, alldata[:, :runtime_box]))
        runtime_pb = alldata[:, :runtime_pb]
        time_plot = plot(xlabel="method", ylabel="runtime")
        boxplot!(time_plot, ["BB" "Box($(box_optimal))" "PB"], [runtime_bb runtime_box runtime_pb], 
            leg=false, 
            linewidth=2,
            colour = [RGB(122/255, 200/255, 255/255) RGB(0/255, 71/255, 119/255) RGB(207/255, 159/255, 205/255) RGB(210/255, 22/255, 53/255)],
            linecolour= :match,
            fillalpha = 0.4)
            savefig(time_plot, "$(filename)_boxplot.pdf")
    end
    if objectives == true && times==true
        combi = plot(obj_plot, time_plot, layout=(1,2), legend=false)
        savefig(combi, "$(filename)_boxplot.pdf")
    end
end