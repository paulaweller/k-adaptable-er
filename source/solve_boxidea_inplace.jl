using JuMP, Gurobi
const GRB_ENV_box_inplace = Gurobi.Env()

"""
    solve_boxes_inplace(K, inst)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_boxes_inplace(K::Int, inst::AllocationInstance)
    time_start = now()
    runtime = 0

    scenario_modell = build_scenario_based_box(inst, K)
    theta, x, y, s, xi = solve_scenario_based_box(scenario_modell, time_start)

    separation_modell = build_separation_problem_box(inst, K, xi)
    zeta, d = solve_separation_problem_box(separation_modell, time_start)
    iteration = 0 # iteration counter

    while zeta > 10^(-6) && (runtime <= 240)
        iteration = iteration + 1

        scenario_modell = update_scenario_based_box!(scenario_modell, K, d)
        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        theta, x, y, s, xi = solve_scenario_based_box(scenario_modell, time_start)

        # find violations
        separation_modell = update_separation_problem_box!(separation_modell, xi)
        zeta, d = solve_separation_problem_box(separation_modell, time_start)
        runtime = (now()-time_start).value/1000
    end
    
    return x, y, s, xi, theta, iteration, runtime

end

"""
    build_scenario_based_boxes(inst, K)

Build the scenario-based K_adaptable problem.
"""
function build_scenario_based_box(inst::AllocationInstance, K::Int)
    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 1)
    J = size(loc_J, 1)
    #c = reshape([norm(loc_I[i,:]-loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    D = inst.D
    pc = inst.pc
    W = inst.W
    

    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_box_inplace); add_bridges = false)
    set_optimizer_attribute(rm, "OutputFlag", 0)
    #set_string_names_on_creation(rm, false) # disable string names for performance improvement

    @expression(rm, c[i=1:I,j=1:J], norm(loc_I[i,:]-loc_J[j,:])); # transportation costs
    @expression(rm, slack_coeff, 10*max(c...))                  # coefficient for slack variables in objective

    @variable(rm, 0 <= w[1:I] <= W, Int)            # first-stage decision
    @variable(rm, 0 <= q[1:I,1:J,1:K] <= W, Int)    # Second stage, wait-and-see decision how to distribute and slack
    @variable(rm, 0 <= s[1:J, 1:K] <= D, Int)       # One set of variables per cell
    @variable(rm, 0 <= ξ[1:J, 1:K] <= D)            # Worst-case scenarios for each box 
    v = rm[:v] = @variable(rm, [1:1, 1:K], binary = true, base_name = "v")                 # which box does scenario t belong to
    
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

function update_scenario_based_box!(scenario_model, K::Int, scenario)

    t = size(scenario_model[:v],1)+1
    J = size(scenario_model[:s],1)
    v_new = @variable(scenario_model, [[t], 1:K], binary = true, base_name = "v")
    scenario_model[:v] = [scenario_model[:v]; v_new]

    @constraint(scenario_model, sum(scenario_model[:v][t,k] for k=1:K) >= 1)        # every demand scenario must be covered by at least one plan
    @constraint(scenario_model, [j=1:J,k=1:K], scenario[j]*scenario_model[:v][t,k] <= scenario_model[:ξ][j,k])   # if plan k covers scenario tau[t], it must be componentwise larger

    return scenario_model
end

function solve_scenario_based_box(scenario_model, time_start)
    # calculate remaining time before cutoff
    time_remaining = 240 + (time_start - now()).value/1000
    # set solver time limit accordingly
    set_time_limit_sec(scenario_model, max(time_remaining,0))
    # solve
    optimize!(scenario_model)
    theta = getvalue(scenario_model[:obj])
    x = round.(Int,getvalue.(scenario_model[:w]))
    y = round.(Int,getvalue.(scenario_model[:q]))
    s = round.(Int,getvalue.(scenario_model[:s]))
    xi = getvalue.(scenario_model[:ξ])
    return theta, x, y, s, xi
end
######################################################################################################################################
function build_separation_problem_box(inst::AllocationInstance, K::Int, ξ_value)

    loc_J = inst.loc_J
    J = size(loc_J, 1)
    D = inst.D
    pc = inst.pc
    us = Model(() -> Gurobi.Optimizer(GRB_ENV_box_inplace); add_bridges = false)
    set_optimizer_attribute(us, "OutputFlag", 0)
    #set_string_names_on_creation(us, false) # disable string names for performance improvement
    
    @variable(us, zeta)     # amount of violation
    @variable(us, 0<= d[1:J] <=D)   # demand scenario // TODO: does it need to be Int?
    @variable(us, z[1:K,1:J], Bin)   # violation indicator
    @variable(us, xi[1:J, 1:K])
    fix.(xi, ξ_value; force = true)

    # d must be in the uncertainty set
    @constraint(us, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))   # bound on aggregated demand
    for (j1,j2) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(us, d[j1]-d[j2] <= norm(loc_J[j1,:]-loc_J[j2,:],Inf))
    end

    @expression(us, M, 2*D+1)
    @constraint(us, [k=1:K], sum(z[k,j] for j in 1:J) == 1)
    #@constraint(us, [k=1:K, j=1:J], zeta + ξ[j,k] <= d[j]*z[k,j])
    @constraint(us, [k=1:K, j=1:J], zeta + M*z[k,j] <= M + d[j] - xi[j,k])

    @objective(us, Max, zeta)
    return us
end

function update_separation_problem_box!(sepmodel, ξ_val)

    fix.(sepmodel[:xi], ξ_val; force = true)
    return sepmodel
end

function solve_separation_problem_box(sepmodel, time_start)    
    # calculate remaining time before cutoff
    time_remaining = 240 + (time_start - now()).value/1000
    # set solver time limit accordingly
    set_time_limit_sec(sepmodel, max(time_remaining,0))

    optimize!(sepmodel)

    return round.(value.(sepmodel[:zeta]), digits = 4), round.(value.(sepmodel[:d]), digits = 2)
end