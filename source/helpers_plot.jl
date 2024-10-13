using LinearAlgebra, Dates, Random, Statistics, JuMP, Plots,StatsPlots
#StatsPlots,


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

    cp = palette(:Set1_3,13)
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
            obj_plot = plot(ylabel="objective", title="BB infeasible", palette=cp)
            if size(θ_box_term,1) > 0
                boxplot!(obj_plot, ["Box t. $(length(θ_box_term))" "Box t. obs"], [θ_box_term obs_box_term], leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4)
                boxplot!(obj_plot, ["PB" "PB obs"], [θ_pb_term obs_pb_term], leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4)
            end
        elseif rel_to_pb == true
            obj_plot = plot(ylabel="relative improvement %", title="BB infeasible",palette=cp)
            if size(θ_box_term,1) > 0
                boxplot!(obj_plot, ["Box t. $(length(θ_box_term))" "Box t. obs"], [(100 .*(θ_pb_term.-θ_box_term)./θ_pb_term) (100 .*(obs_pb_term.-obs_box_term)./obs_pb_term)], leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4)
            end
        end
        if size(θ_box_unterm,1) > 0
            if rel_to_pb == false
                boxplot!(obj_plot, ["Box unt. $(length(θ_box_unterm))" "Box unt. obs"], [θ_box_unterm obs_box_unterm], leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4)
                boxplot!(obj_plot, ["PB" "PB obs"], [θ_pb_unterm obs_pb_unterm], leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4)
            elseif rel_to_pb == true
                boxplot!(obj_plot, ["Box unt. $(length(θ_box_unterm))" "Box unt. obs"], [(100 .*(θ_pb_unterm.-θ_box_unterm)./θ_pb_unterm) (100 .*(obs_pb_unterm.-obs_box_unterm)./obs_pb_unterm)], 
                leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4)
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
                markerstrokewidth=0,
                fillalpha = 0.4)
        end
        if length(time_box_unterm) > 0
            boxplot!(time_plot, ["Box unt." "PB"], [time_box_unterm time_pb_unterm], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                markerstrokewidth=0,
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
                boxplot!(obj_plot_bb_box, ["BB" "BB obs" "Box" "Box obs" "PB" "PB obs"], [θ_bb obs_bb θ_box obs_box θ_pb obs_pb], 
                colour=[cp[1] cp[3] cp[6] cp[8] cp[11] cp[13]],leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4) #TODO relative
                savefig(obj_plot_bb_box, "source/plots/termination/k$(K)_$(terminated)_terminate_obj.pdf")
            elseif rel_to_pb == true
                obj_plot_bb_box = plot(ylabel="relative improvement %", title="$(nrow(plot_data)) instances", palette=cp)
                boxplot!(obj_plot_bb_box, ["BB" "BB obs" "Box" "Box obs"], [(100 .*(θ_pb.-θ_bb)./θ_pb) (100 .*(obs_pb.-obs_bb)./obs_pb) (100 .*(θ_pb.-θ_box)./θ_pb) (100 .*(obs_pb.-obs_box)./obs_pb)], 
                leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4) #TODO relative
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
                markerstrokewidth=0,
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
    cp = palette(:Set1_3,15)
    if method == "bb"
        color5 = [cp[1] cp[2] cp[3] cp[4] cp[5]]
        color3 = [cp[1] cp[2] cp[3]]
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
        color5 = [cp[6] cp[7] cp[8] cp[9] cp[10]]
        color3 = [cp[6] cp[7] cp[8]]
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
        color5 = [cp[11] cp[12] cp[13] cp[14] cp[15]]
        color3 = [cp[13] cp[14] cp[15]]
        term = alldata_dirty
        observable ? (θ_key = :obs_pb) : (θ_key = :θ_pb)
        time_key = :runtime_pb
    else
        error("method is invalid, must be one of bb, box, pb")
    end
    # sort by k
    term_k1 = term[term[!,:k] .== 1,:]
    term_k2 = term[term[!,:k] .== 2,:]
    term_k3 = term[term[!,:k] .== 3,:]

    if observable 
        term_is1 = term_k1
        term_is2 = term_k2
        term_is3 = term_k3
        term_is4 = term[term[!,:k] .== 4,:]
        term_is5 = term[term[!,:k] .== 5,:]
        sort!(term_is1, [:n, :m, :pc, :instance])
        sort!(term_is2, [:n, :m, :pc, :instance])
        sort!(term_is3, [:n, :m, :pc, :instance])
        sort!(term_is4, [:n, :m, :pc, :instance])
        sort!(term_is5, [:n, :m, :pc, :instance])
    else
        term_is1 = semijoin(term_k1, term_k2, on = [:instance, :n, :m, :pc])
        term_is1 = semijoin(term_is1, term_k3, on = [:instance, :n, :m, :pc])
        sort!(term_is1, [:n, :m, :pc, :instance])

        term_is2 = semijoin(term_k2, term_k1, on = [:instance, :n, :m, :pc])
        term_is2 = semijoin(term_is2, term_k3, on = [:instance, :n, :m, :pc])
        sort!(term_is2, [:n, :m, :pc, :instance])

        term_is3 = semijoin(term_k3, term_k1, on = [:instance, :n, :m, :pc])
        term_is3 = semijoin(term_is3, term_k2, on = [:instance, :n, :m, :pc])
        sort!(term_is3, [:n, :m, :pc, :instance])
    end
    if relative == true 
        obj_plot = plot(ylabel="relative improvement of objective",palette=cp)
        θ1 = 100 .*(1 .- term_is1[:, θ_key]./term_is1[:, θ_key])
        θ2 = 100 .*(1 .- term_is2[:, θ_key]./term_is1[:, θ_key])
        θ3 = 100 .*(1 .- term_is3[:, θ_key]./term_is1[:, θ_key])
        if observable == true
            θ4 = 100 .*(1 .- term_is4[:, θ_key]./term_is1[:, θ_key])
            θ5 = 100 .*(1 .- term_is5[:, θ_key]./term_is1[:, θ_key])
        end
    elseif relative == false
        obj_plot = plot(ylabel="objective",palette=cp)
        θ1 = term_is1[:, θ_key]
        θ2 = term_is2[:, θ_key]
        θ3 = term_is3[:, θ_key]
        if observable == true
            θ4 = term_is4[:, θ_key]
            θ5 = term_is5[:, θ_key]
        end
    else 
        error("value for relative is invalid")
    end
    n_inst = length(θ2)

    if observable
        boxplot!(obj_plot, ["K=1" "K=2" "K=3" "K=4" "K=5"], [θ1 θ2 θ3 θ4 θ5], 
        color=color5,
        leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,markeropacity=1,fillalpha = 0.4, ylims=(-15,60))
    else
        boxplot!(obj_plot, ["K=1" "K=2" "K=3"], [θ1 θ2 θ3], 
        color=color3,leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4)
    end
    if relative == true 
        observable ? savefig(obj_plot, "source/plots/k_comparison/k_evol_obsrelobj_k_$(method)_$(n_inst).pdf") : savefig(obj_plot, "source/plots/k_comparison/relobj_k_$(method)_$(n_inst).pdf")
    else 
        observable ? savefig(obj_plot, "source/plots/k_comparison/obsobj_k_$(method)_$(n_inst).pdf") : savefig(obj_plot, "source/plots/k_comparison/obj_k_$(method)_$(n_inst).pdf")
    end
    if observable == false 
        t1 = alldata[alldata[!,:k] .== 1, time_key]##term_is1[:, time_key]
        t2 = alldata[alldata[!,:k] .== 2, time_key]##term_is2[:, time_key]
        t3 = alldata[alldata[!,:k] .== 3, time_key]##term_is3[:, time_key]
        time_plot = plot(ylabel="runtime",palette=cp)
        violin!(time_plot, ["K=1" "K=2" "K=3"], [t1 t2 t3], 
            leg=false, 
            xtickfontsize=12,
            ytickfontsize=12,
            xguidefontsize=14,
            yguidefontsize=14,
            linewidth=2,
            linecolour= :match,
            markerstrokewidth=0,
            color=color3,
            fillalpha = 0.4,palette=cp)

        savefig(time_plot, "source/plots/k_comparison/violin_time_k_$(method)_all.pdf")
        # xtickfontsize=18,
        # ytickfontsize=18,
        # xguidefontsize=18,
        # yguidefontsize=18,
        # legendfontsize=18
    end
    # combi = plot(obj_plot_bb_box, time_plot, layout=(1,2), legend=false)
    # savefig(combi, "$(filename)_boxplot.pdf")

end

function plot_pc_vs_time(filename; time=false, objective=false, K=[1,2,3,4,5])
    alldata_all = DataFrame(CSV.File("$(filename).csv"))
    alldata = alldata_all[in.(alldata_all[!,:k],(K,)),:]
    cp = palette(:Set1_6)
    # replace!(alldata.runtime_bb, "i" => "3600.0")
    # alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
    
    data1 = alldata[alldata[!,:pc] .== 0.1,:]
    data3 = alldata[alldata[!,:pc] .== 0.3,:]

    if time == true
        time_plot1 = plot(xlabel="pc = 0.1",palette=cp)
        boxplot!(time_plot1, ["bb" "box" "pb"], [data1[!,:runtime_bb] data1[!,:runtime_box] data1[!,:runtime_pb]], 
            leg=false, 
            linewidth=2,
            linecolour= :match,
            markerstrokewidth=0,
            fillalpha = 0.4)

            time_plot3 = plot(xlabel="pc = 0.3",palette=cp)
            boxplot!(time_plot3, ["bb" "box" "pb"], [data3[!,:runtime_bb] data3[!,:runtime_box] data3[!,:runtime_pb]], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                markerstrokewidth=0,
                fillalpha = 0.4)

        cplot = plot(time_plot1, time_plot3, ylabel="runtime")
        savefig(cplot, "source/plots/percentage/pc_time_$(K).pdf")
    end
    if objective == true
        alldata_clean = alldata[alldata[!,:θ_bb] .!= 1e20,:]
        # alldata_clean = alldata[alldata[!,:it_bb] .!= "s",:]
        # alldata_clean[!,:θ_bb] = parse.(Float64, alldata_clean[!,:θ_bb])

        data1 = alldata_clean[alldata_clean[!,:pc] .== 0.1,:]
        data3 = alldata_clean[alldata_clean[!,:pc] .== 0.3,:]

        data1 = semijoin(data1, data3, on = [:instance, :n, :m])
        data3 = semijoin(data3, data1, on = [:instance, :n, :m])

        time_plot1 = plot(xlabel="pc = 0.1",palette=cp)
        boxplot!(time_plot1, ["bb" "box" "pb"], [data1[!,:θ_bb] data1[!,:θ_box] data1[!,:θ_pb]], 
            leg=false, 
            linewidth=2,
            linecolour= :match,
            markerstrokewidth=0,
            fillalpha = 0.4)

            time_plot3 = plot(xlabel="pc = 0.3",palette=cp)
            boxplot!(time_plot3, ["bb" "box" "pb"], [data3[!,:θ_bb] data3[!,:θ_box] data3[!,:θ_pb]], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                markerstrokewidth=0,
                fillalpha = 0.4)

        cplot = plot(time_plot1, time_plot3, ylabel="objective")
        savefig(cplot, "source/plots/percentage/pc_objective.pdf")
    end
end

function plot_evol(zetas, filename; xlimits=[0,3600], relative=true, name="zeta", last=false)
    cp = palette(:Set1_6)
    if name[1:2] == "bb"
        mycolor = cp[1]
    elseif name[1:2] == "bo"
        mycolor = cp[2]
    else
        mycolor = cp[3]
    end
    time_plot = plot(xlims=xlimits,palette=cp)
    # extract pb data for comparison
    df = DataFrame(CSV.File("$(filename).csv"))
    for (key, values) in zetas
        # key is [n,m,pc,k,l]
        pb = df[ ( df.n .== key[1] ) .& ( df.m .== key[2] ) .& (df.pc .== key[3]) .& (df.k .== key[4]) .& (df.instance .== key[5]), :θ_pb]

        # Extract x and y coordinates from each 2D vector
        if length(values) > 1 && values != "Vector{Float64}[]"
            coordinates = vecfromstr(values)
            
            # Extract x and y coordinates from each 2D vector
            x_vals = coordinates[:,1]
            if name != "zeta"
                relative ? (y_vals = (pb .- coordinates[:,2])./pb) : (y_vals = coordinates[:,2])
            else
                if coordinates[1,2] == 0
                    if last == true
                        relative ? (y_vals = (coordinates[end,2] .- coordinates[:,2])./coordinates[end,2]) : (y_vals = coordinates[:,2])
                    else 
                        relative ? (y_vals = (coordinates[2,2] .- coordinates[:,2])./coordinates[2,2]) : (y_vals = coordinates[:,2])
                    end
                else 
                    relative ? (y_vals = (coordinates[1,2] .- coordinates[:,2])./coordinates[1,2]) : (y_vals = coordinates[:,2])
                end
            end
            # Plot the data
            plot!(time_plot, x_vals[1:(end)], y_vals[1:(end)], 
            label="",line = (0.1,1), color=mycolor, ylims=(0,0.4),xlabel="seconds",ylabel="objective improvement relative to PB")
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
    cp = palette(:Set1_6)
    ms = 5
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
            
            scatter!(obj_plot, (36/(n_interval)).*lastk, z[lastk], label="k=$k", markerstrokewidth=0, markersize=2+1.5*ms, markeropacity=0.6, palette=cp)
            #(i==2) ? scatter!(obj_plot, (36/(n_interval)).*lastk, z[lastk], label="", markerstrokewidth=0, markersize=10, markeralpha=0.3, palette=cp) : nothing
            ms = ms-1
        end
    if status =="terminated"
        relative ? savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_termin_rel_$(K).pdf") : savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_termin_$(K).pdf")
    else 
        relative ? savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_unterm_rel_$(K).pdf") : savefig(obj_plot, "source/plots/zetas/last_zeta/zeta_dist_untermin_$(K).pdf")
    end
end

number_in_interval(lb, rb, zetas) = count(i-> ((i <= rb) && (i > lb)), zetas)


function plot_size_vs_time(filename; K=[2])
    alldata = DataFrame(CSV.File("$(filename).csv"))

    # replace!(alldata.runtime_bb, "i" => "3600.0")
    # alldata[!,:runtime_bb] = parse.(Float64, alldata[!,:runtime_bb])
    cp = palette(:Set1_6)
    data = alldata#[alldata[!,:pc] .== percentage,:]
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
            violin!(time_plot, ["BB" "Box"], [data_k[:,:runtime_bb] data_k[:,:runtime_box]], 
                leg=false, linewidth=2,linecolour= :match,markerstrokewidth=0,fillalpha = 0.4,palette=cp)
            push!(plots, time_plot)
        end
    end
    cplot = plot(plots..., layout=(3,3))
    savefig(cplot, "source/plots/size/violin_runtime_pcall_K$(K).pdf")
end

function observable_plot(o_pb, o_bb, o_box, oo_pb, oo_box,k, pc; rel=true)
    cp = palette(:Set1_6)
    if rel == true
        plot1 = plot(ylabel="observable improvement",palette=cp)
        boxplot!(plot1, ["bb ($(length(o_bb)))" "box"], [100 .*(o_pb.-o_bb)./o_pb 100 .*(o_pb.-o_box)./o_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                markerstrokewidth=0,
                fillalpha = 0.4)
        if length(oo_pb) > 0
            boxplot!(plot1, ["bb infeas ($(length(oo_box)))"], [100 .*(oo_pb.-oo_box)./oo_pb], 
                    leg=false, 
                    linewidth=2,
                    linecolour= :match,
                    markerstrokewidth=0,
                    fillalpha = 0.4)
        end
        savefig(plot1, "source/plots/observable/rel_objective_$(pc)_k$(k).pdf")
    elseif rel == false
        plot1 = plot(ylabel="observable worst-case objective",palette=cp)
        boxplot!(plot1, ["bb($(length(o_bb)))" "box" "pb"], [o_bb o_box o_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                markerstrokewidth=0,
                fillalpha = 0.4)
        if length(oo_pb) > 0
                boxplot!(plot1, ["bb infeas ($(length(oo_box))): box" " bb infeas: pb"], [oo_box oo_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                markerstrokewidth=0,
                fillalpha = 0.4)
        end
        savefig(plot1, "source/plots/observable/objective_$(pc)_k$(k).pdf")
    else error("optinal argument rel must be true or false")
    end
end
