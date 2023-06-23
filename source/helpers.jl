using LinearAlgebra, Dates, Random, Statistics, JuMP
#StatsPlots,
struct AllocationInstance
    loc_I::Array{Float64,2}
    loc_J::Array{Float64,2}
    W::Float64
    D::Float64
    pc::Float64
end
#const GRB_ENV = Gurobi.Env()
"""
    generate_instance(I_inst, J_inst, seed, demand_bound=5, cont_perc=0.5, agg_supply_bound=round(Int, cont_perc*demand_bound*J), plot_loc=false)

Generate a problem instance with the given parameters. Returns: loc_I_inst, loc_J_inst, demand_bound, cont_perc, agg_supply_bound

Fields:

    I_inst              Number of service points
    J_inst              Number of demand points
    seed                For reproducibility
    demand_bound        Upper bound of demand at one demand point
    cont_perc           Percentage of damage caused by the contingency
    agg_supply_bound    Aggregated supply bound, default: maximal aggregated demand 
    plot_loc            Should the locations of service and demand points be plotted?
    loc_max             How large is the grid
"""
function generate_instance(I_inst, J_inst, seed; demand_bound=5, cont_perc=0.5, agg_supply_bound=round(Int, cont_perc*demand_bound*J_inst), plot_loc=false, loc_max=10)

    # seed for reproducibility
    Random.seed!(seed)
    
    # set of locations for service and demand points
    loc_set = Vector{Float64}[]
    # locations of service points
    loc_I_inst = Vector{Float64}[]
    # locations of demand points
    loc_J_inst = Vector{Float64}[]

    # no two points can have the same location
    while length(loc_set) !== I_inst+J_inst  
        
        # randomly generate locations of service points as x,y-coordinates
        loc_I_inst = rand(1:loc_max,2,I_inst)
        # randomly generate locations of demand points as coordinates
        loc_J_inst = rand(1:loc_max,2,J_inst)
        # add to set (not list) of locations, i.e. a location won't be listed twice
        loc_set = union(union(loc_I_inst[:,i] for i in 1:I_inst),union(loc_J_inst[:,j] for j in 1:J_inst))
    end

    # sort locations by x-coordinate for simplicity
    loc_I_inst = sort(loc_I_inst, dims=2)
    loc_J_inst = sort(loc_J_inst, dims=2)

    instance = AllocationInstance(loc_I_inst,loc_J_inst,agg_supply_bound,demand_bound,cont_perc)
    # plot service and demand points
    if plot_loc == true

        scatter(loc_I_inst[1,:],loc_I_inst[2,:],
                    lims=[0,loc_max+1.2],
                    series_annotations = text.(1:J_inst, :inside),
                    markersize = 20,
                    lab="service point")
        display(scatter!(loc_J_inst[1,:],loc_J_inst[2,:], 
                            series_annotations = text.(1:J_inst, :inside),
                            shape = :square,
                            markersize = 15,
                            lab="demand point",
                            aspect_ratio=:equal
                            ))
    end

    return instance
end

function set_remaining_time(model::Model, time_start::DateTime)
    # calculate remaining time before cutoff
    time_remaining = 240 + (time_start - now()).value/1000
    # set solver time limit accordingly
    set_time_limit_sec(model, max(time_remaining,0))
end

function branch_partition!(N::Vector{Vector{Vector{Vector{Float64}}}}, tau::Vector{Vector{Vector{Float64}}}, xi::Array{Float64,1}, knew::Int64)
    for k in 1:knew
        # each child node is the current uncset configuration...
        tau_temp = copy(tau)
        # ...with the new scenario added to uncset number k
        # tau_temp[k] = union(tau_temp[k], [xi])
        unique!(push!(tau_temp[k], xi))
        #N = union(N, [tau_temp])
        push!(N, tau_temp)  
    end
    nothing #return N
end

function number_of_childnodes(tau::Vector{Vector{Vector{Float64}}})
    # sort the partition by size of the subsets
    sort!(tau, by= x-> size(x), rev = true)
    #define Knew the first empty subset
    Knew = length(tau)
    if isempty(tau[end]) 
        Knew = findfirst(x-> isempty(x), tau)
    end
    return Knew
end

"""
    find_plan(d_real, q, c)

Find the best of the plans in q for this uncertainty realization d_real.
"""
function find_plan(d_real, q, c)
    c_slack = 1+max(c...)
    obj_val = 10^10
    best_q = []
    for p in 1:size(q, 3)
        obj_val_it = c_slack*sum(max(0, d_real[j]-sum(q[i,j,p] for i in 1:size(q,1))) for j in 1:size(q,2)) + sum(c[i,j]*q[i,j,p] for i in 1:size(q,1), j in 1:size(q,2))
        if obj_val_it < obj_val
            obj_val = obj_val_it
            best_q = q[:,:,p]
        end
    end
    return obj_val, best_q
end

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

function box_plot_from_files(filenames, labelnames)

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