using JuMP, Gurobi
const GRB_ENV_box = Gurobi.Env()

"""
    solve_boxes(K, instance)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_boxes(K::Int64, inst::AllocationInstance; time_limit::Float64 = 240.0)
    time_start = now()
    runtime = 0.0
    tau = Vector{Float64}[]            # set of scenarios 

    theta, x, y, s, xi = solve_scenario_based_boxes(tau, inst, K, time_start, time_limit)
    zeta, d = solve_separation_problem_boxes(inst, K, xi, time_start, time_limit)
    iteration = 0 # iteration counter

    
    while zeta > 1e-6 && (runtime <= 120.0)
        iteration = iteration + 1
        push!(tau, ceil.(round.(d, digits=4)))
        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        theta, x, y, s, xi = solve_scenario_based_boxes(tau, inst, K, time_start, time_limit)

        # find violations
        zeta, d = solve_separation_problem_boxes(inst, K, xi, time_start, time_limit)
        runtime = (now()-time_start).value/1000
    end
    
    return x, y, s, xi, theta, iteration, runtime

end

"""
    solve_scenario_based_boxes(tau, inst, K, starttime)

Solve the scenario-based K_adaptable problem for the uncertainty sets tau.
"""
function solve_scenario_based_boxes(tau::Vector{Vector{Float64}}, inst::AllocationInstance, K::Int64, time_start::DateTime, time_limit::Float64)
    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 2)
    J = size(loc_J, 2)
    D = inst.D
    W = inst.W
    T = length(tau)

    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_box); add_bridges = false)
    set_optimizer_attribute(rm, "OutputFlag", 0)
    set_string_names_on_creation(rm, false) # disable string names for performance improvement
    set_optimizer_attribute(rm, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)

    @expression(rm, c[i=1:I,j=1:J], norm(loc_I[:,i]-loc_J[:,j])); # transportation costs
    @expression(rm, slack_coeff, 10*max(c...))                  # coefficient for slack variables in objective


    @variable(rm, 0 <= w[1:I] <= W, Int)            # first-stage decision
    @variable(rm, 0 <= q[1:I,1:J,1:K] <= W, Int)    # Second stage, wait-and-see decision how to distribute and slack
    @variable(rm, 0 <= s[1:J, 1:K] <= D, Int)       # One set of variables per cell
    @variable(rm, 0 <= ξ[1:J, 1:K] <= D)            # Worst-case scenarios for each box 
    @variable(rm, v[1:T, 1:K], Bin)                 # which box does scenario t belong to
    
    @constraint(rm, sum(w[i] for i in 1:I) <= W)                    # supply limit
    @constraint(rm, [i=1:I,k=1:K], sum(q[i,j,k] for j in 1:J) <= w[i])    # service point limit
        
    @variable(rm, 0 <= obj <= 10^10)                        # The objective function will be the maximum of the objective function across all the cells
    @variable(rm, 0<= z[1:K]<=10^10)                        # A variable to track each cells objective value
    @objective(rm, Min, obj)
        
    # Constrain objective function for this cell
    @constraint(rm, [k=1:K], z[k] >= slack_coeff*sum(s[j,k] for j in 1:J) + sum(c[i,j]*q[i,j,k] for i in 1:I, j in 1:J))
    @constraint(rm, [k=1:K], obj >= z[k])

    # Demand 
    @constraint(rm, [j=1:J, k=1:K], sum(q[i,j,k] for i in 1:I)+s[j,k] >= ξ[j,k])


    @constraint(rm, [t=1:T], sum(v[t,k] for k=1:K) >= 1)        # every demand scenario must be covered by at least one plan
    @constraint(rm, [j=1:J,t=1:T,k=1:K], tau[t][j]*v[t,k] <= ξ[j,k])   # if plan k covers scenario tau[t], it must be componentwise larger

    set_remaining_time(rm, time_start, time_limit)
    # solve
    optimize!(rm)
    if result_count(rm) == 0
        return 1e10, zeros(Float64, I), zeros(Float64, I, J, K), zeros(Float64, J, K), zeros(Float64, J,K)
    end
    theta = value(obj)
    x = value.(w)
    y = value.(q)
    s_val = value.(s)
    xi = value.(ξ)
    return theta, x, y, s_val, xi
end

function solve_separation_problem_boxes(inst::AllocationInstance, K::Int64, ξ::Array{Float64,2}, time_start::DateTime, time_limit::Float64)

    loc_J = inst.loc_J
    J = size(loc_J, 2)
    D = inst.D
    pc = inst.pc
    us = Model(() -> Gurobi.Optimizer(GRB_ENV_box); add_bridges = false)
    set_optimizer_attribute(us, "OutputFlag", 0)
    set_string_names_on_creation(us, false) # disable string names for performance improvement
    set_optimizer_attribute(us, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)
    

    @variable(us, zeta)     # amount of violation
    @variable(us, 0<= d[1:J] <=D)   # demand scenario // TODO: does it need to be Int?
    @variable(us, z[1:K,1:J], Bin)   # violation indicator

    # d must be in the uncertainty set
    @constraint(us, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))   # bound on aggregated demand
    for (j2,j1) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(us, d[j1]-d[j2] <= norm(loc_J[:,j1]-loc_J[:,j2],Inf))
    end

    @expression(us, M, 2*D+1)
    @constraint(us, [k=1:K], sum(z[k,j] for j in 1:J) == 1)
    #@constraint(us, [k=1:K, j=1:J], zeta + ξ[j,k] <= d[j]*z[k,j])
    @constraint(us, [k=1:K, j=1:J], zeta + M*z[k,j] <= M + d[j] - ξ[j,k])

    @objective(us, Max, zeta)

    set_remaining_time(us, time_start, time_limit)
    optimize!(us)
    if result_count(us) == 0
        return 0.0, zeros(Float64, J)
    end
    return value.(zeta), value.(d)
end