using JuMP, Gurobi
const GRB_ENV_box = Gurobi.Env()

"""
    solve_box(K, inst)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_box(K::Int64, inst::AllocationInstance; time_limit::Float64 = 240.0)
    time_start = now()
    runtime = 0.0

    scenario_modell = build_scenario_based_box(inst, K)
    theta, x, y, s, xi = solve_scenario_based_box(scenario_modell, time_start, time_limit)

    separation_modell = build_separation_problem_box(inst, K, xi)
    zeta, d = solve_separation_problem_box(separation_modell, time_start, time_limit)
    iteration = 0 # iteration counter

    while zeta > 1e-6 && (runtime <= time_limit)
        iteration = iteration + 1

        update_scenario_based_box!(scenario_modell, K, ceil.(round.(d, digits = 4)))
        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        theta, x, y, s, xi = solve_scenario_based_box(scenario_modell, time_start, time_limit)

        # find violations
        update_separation_problem_box!(separation_modell, xi)
        zeta, d = solve_separation_problem_box(separation_modell, time_start, time_limit)
        runtime = (now()-time_start).value/1000
    end
    
    return x, y, s, xi, theta, iteration, runtime

end

"""
    build_scenario_based_box(inst, K)

Build the scenario-based K_adaptable problem.
"""
function build_scenario_based_box(inst::AllocationInstance, K::Int64)
    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 2)
    J = size(loc_J, 2)
    D = inst.D
    D_max = max(D...)
    W = inst.W
    
    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_box); add_bridges = false)
    set_optimizer_attribute(rm, "OutputFlag", 0)
    set_string_names_on_creation(rm, false) # disable string names for performance improvement
    set_optimizer_attribute(rm, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)
    set_optimizer_attribute(rm, "Threads", 8)

    @expression(rm, c[i=1:I,j=1:J], norm(loc_I[:,i]-loc_J[:,j])); # transportation costs
    @expression(rm, slack_coeff, 1000.0*norm(c,Inf))                  # coefficient for slack variables in objective

    @variable(rm, 0 <= w[1:I] <= W, Int)            # first-stage decision
    @variable(rm, 0 <= q[1:I,1:J,1:K] <= W, Int)    # Second stage, wait-and-see decision how to distribute and slack
    @variable(rm, 0 <= s[1:J, 1:K] <= D_max, Int)       # One set of variables per cell
    @variable(rm, 0 <= ξ[1:J, 1:K] <= D_max)            # Worst-case scenarios for each box 
    # v = rm[:v] = @variable(rm, [1:1, 1:K], binary = true, base_name = "v")                 # which box does scenario t belong to
    
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

    return rm
end

function update_scenario_based_box!(scenario_model::Model, K::Int64, scenario::Vector{Float64})

    # t = size(scenario_model[:v],1)+1
    J = size(scenario_model[:s],1)
    # v_new = @variable(scenario_model, [[t], 1:K], binary = true, base_name = "v")
    # scenario_model[:v] = [sckenario_model[:v]; v_new]

    # @constraint(scenario_model, sum(scenario_model[:v][t,k] for k=1:K) >= 1)        # every demand scenario must be covered by at least one plan
    # @constraint(scenario_model, [j=1:J,k=1:K], scenario[j]*scenario_model[:v][t,k] <= scenario_model[:ξ][j,k])   # if plan k covers scenario tau[t], it must be componentwise larger

    v_new = @variable(scenario_model, [1:K], binary = true)
    # scenario_model[:v] = [sckenario_model[:v]; v_new]

    @constraint(scenario_model, sum(v_new[k] for k=1:K) >= 1)        # every demand scenario must be covered by at least one plan
    @constraint(scenario_model, [k=1:K,j=1:J], scenario[j]*v_new[k] <= scenario_model[:ξ][j,k])   # if plan k covers scenario tau[t], it must be componentwise larger

    nothing #return scenario_model
end

function solve_scenario_based_box(scenario_model::Model, time_start::DateTime, time_limit::Float64)
    set_remaining_time(scenario_model, time_start, time_limit)
    optimize!(scenario_model)
    if result_count(scenario_model) == 0
        I, J, K = size(scenario_model[:q])
        return 1e10, zeros(Float64, I), zeros(Float64, I, J, K), zeros(Float64, J, K), zeros(Float64, J,K)
    end
    theta::Float64 = value(scenario_model[:obj])
    x::Vector{Float64} = value.(scenario_model[:w])
    y::Array{Float64,3} = value.(scenario_model[:q])
    s_val::Array{Float64,2} = value.(scenario_model[:s])
    xi::Array{Float64,2} = value.(scenario_model[:ξ])
    return theta, x, y, s_val, xi
end
######################################################################################################################################
function build_separation_problem_box(inst::AllocationInstance, K::Int64, ξ_value::Array{Float64,2})

    loc_J = inst.loc_J
    J = size(loc_J, 2)
    D = inst.D
    D_max = max(D...)
    pc = inst.pc
    us = Model(() -> Gurobi.Optimizer(GRB_ENV_box); add_bridges = false)
    set_optimizer_attribute(us, "OutputFlag", 0)
    set_string_names_on_creation(us, false) # disable string names for performance improvement
    set_optimizer_attribute(us, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)
    set_optimizer_attribute(us, "Threads", 8)
    
    @variable(us, zeta)     # amount of violation
    @variable(us, 0<= d[1:J])
    @constraint(us, [j=1:J], d[j] <= D[j])   # demand scenario // TODO: does it need to be Int?
    @variable(us, z[1:K,1:J], Bin)   # violation indicator
    @variable(us, xi[1:J, 1:K])
    fix.(xi, ξ_value; force = true)

    # d must be in the uncertainty set
    @constraint(us, sum(d[j] for j in 1:J) <= floor(pc*sum(D)))   # bound on aggregated demand
    for (j2,j1) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(us, d[j1]-d[j2] <= norm(loc_J[:,j1]-loc_J[:,j2],Inf))
    end

    @expression(us, M, 2*D_max+1)
    @constraint(us, [k=1:K], sum(z[k,j] for j in 1:J) == 1)
    #@constraint(us, [k=1:K, j=1:J], zeta + ξ[j,k] <= d[j]*z[k,j])
    @constraint(us, [k=1:K, j=1:J], zeta + M*z[k,j] <= M + d[j] - xi[j,k])

    @objective(us, Max, zeta)
    return us
end

function update_separation_problem_box!(sepmodel::Model, ξ_val::Array{Float64,2})
    fix.(sepmodel[:xi], ξ_val; force = true)
    nothing #return sepmodel
end

function solve_separation_problem_box(sepmodel::Model, time_start::DateTime, time_limit::Float64)    
    set_remaining_time(sepmodel, time_start, time_limit)
    optimize!(sepmodel)
    if result_count(sepmodel) == 0
        J = length(sepmodel[:d])
        return 0.0, zeros(Float64, J)
    end
    zeta_val::Float64 = value(sepmodel[:zeta])
    d_val::Vector{Float64} = value.(sepmodel[:d])
    return zeta_val, d_val
end