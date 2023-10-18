using LinearAlgebra, Dates, Random, Statistics, JuMP

include("solve_iter.jl")
include("solve_bb.jl")
include("solve_boxes_inplace.jl")
include("solve_bb_inplace.jl")
include("solve_boxes.jl")

"""
    run_instance(k, problem_instance; tlim=250.0, pb=true, box=false, box_inplace=true, bb=true, bb_inplace=false)

Calculate k-adaptabe solution for an instance of the preallocation problem, using the specified algorithms. Returns results as a dictionary.
"""
function run_instance(k, problem_instance; tlim=250.0, pb=true, box=false, box_inplace=true, bb=true, bb_inplace=false)

    results = Dict()  # for storing single-value results 

    println("\n")
    if pb==true
        # how many iterations are necessary for a k-adaptable solution?
        it = round(Int, 1 + (k-1)/(dn-1))
        # run partition-and-bound method
        println("Starting partition-and-bound...")
        theta_pb, x_pb, y_pb, p_pb, p_true_pb, runtime_pb = k_adapt_solution(it, problem_instance)  # TODO rename method to partition-and-bound AND time limit?
        results_pb = Dict(
            :θ_pb       => theta_pb, 
            # :x_pb       => x_pb, 
            # :y_pb       => y_pb, 
            # :p_pb       => p_pb, 
             :p_true_pb  => p_true_pb, 
            :runtime_pb => runtime_pb
        )
        merge!(results, results_pb)
        println("Finished partition-and-bound. Runtime: ", runtime_pb, "s\n")
    end
    if box==true
        println("Starting box-and-cut...")
        x_box, y_box, s_box, xi_box, theta_box, it_box, runtime_box = solve_boxes(k, problem_instance) # TODO time limit?
        results_box = Dict(
            # :x_box       => x_box, 
            # :y_box       => y_box, 
            # :s_box       => s_box, 
            # :ξ_box       => xi_box, 
            :θ_box       => theta_box, 
            :it_box      => it_box, 
            :runtime_box => runtime_box
        )
        merge!(results, results_box)
        println("Finished box-and-cut. Runtime: ", runtime_box, "s\n")
    end
    if box_inplace==true
        println("Starting box-and-cut inplace version...")
        x_box, y_box, s_box, xi_box, theta_box, it_box, runtime_box = solve_boxes_inplace(k, problem_instance, time_limit = tlim)
        results_box_inplace = Dict(
            # :x_box_inplace       => x_box, 
            # :y_box_inplace       => y_box, 
            # :s_box_inplace       => s_box, 
            # :ξ_box_inplace       => xi_box, 
            :θ_box_inplace       => theta_box, 
            :it_box_inplace      => it_box, 
            :runtime_box_inplace => runtime_box
        )
        merge!(results, results_box_inplace)
        println("Finished box-and-cut inplace version. Runtime: ", runtime_box, "s\n")
    end
    if bb==true
        println("Starting branch-and-bound...")
        x_general, y_general, s_general, partition, theta_general, it_general, runtime_general = solve_bb_general(k, problem_instance, time_limit = tlim) # TODO save partition?
        results_bb = Dict(
            # :x_bb           => x_general, 
            # :y_bb           => y_general, 
            # :s_bb           => s_general, 
            # :partition_bb   => partition,
            :θ_bb           => theta_general, 
            :it_bb          => it_general, 
            :runtime_bb     => runtime_general
        )
        merge!(results, results_bb)
        println("Finished branch-and-bound. Runtime: ", runtime_general, "s\n")
    end
    if bb_inplace==true
        println("Starting branch-and-bound inplace version...")
        x_gen2, y_gen2, s_gen2, partition_gen2, theta_gen2, it_gen2, runtime_gen2 = solve_bb_inplace(k, problem_instance, time_limit = tlim) # TODO save partition?
        results_bb_inplace = Dict(
            # :x_bb_inplace           => x_gen2, 
            # :y_bb_inplace           => y_gen2, 
            # :s_bb_inplace           => s_gen2,
            # :partition_bb_inplace   => partition_gen2, 
            :θ_bb_inplace           => theta_gen2, 
            :it_bb_inplace          => it_gen2, 
            :runtime_bb_inplace     => runtime_gen2
        )
        merge!(results, results_bb_inplace)
        println("Finished branch-and-bound inplace version. Runtime: ", runtime_gen2, "s\n")
    end
    return results
end

function set_remaining_time(model::Model, time_start::DateTime, time_limit::Float64)
    # calculate remaining time before cutoff
    time_remaining = time_limit + (time_start - now()).value/1000
    # set solver time limit accordingly
    set_time_limit_sec(model, max(time_remaining,0))
end

function branch_partition!(N::Vector{Vector{Vector{Vector{Float64}}}}, tau::Vector{Vector{Vector{Float64}}}, xi::Array{Float64,1}, knew::Int64)
    for k in 1:knew
        # each child node is the current uncset configuration...
        tau_temp = copy(tau)
        # ...with the new scenario added to uncset number k
        tau_temp_k = copy(tau_temp[k])
        # tau_temp[k] = union(tau_temp[k], [xi])
        tau_temp[k] = vcat(tau_temp_k, [xi])
        unique!(tau_temp[k])
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