using LinearAlgebra, Dates, Random, Statistics, JuMP

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