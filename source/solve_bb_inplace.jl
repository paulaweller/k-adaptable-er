using JuMP, Gurobi
const GRB_ENV_bb_inplace = Gurobi.Env()

"""
    solve_bb_inplace(K, inst)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_bb_inplace(K::Int64, inst::AllocationInstance)
    time_start = now()  # measure start of runtime
    runtime = 0.0         # initiate variable for runtime
    # set of unexplored nodes, each node consists of K (disjoint) subsets of the uncertainty set
    # we start with K empty sets
    N = [[Iterators.repeated(Vector{Int64}[],K)...]]    # N[a,b,c] = N[partition,subset,demandpoint]
    # the incumbent
    theta_i = 1e10     # objective value
    x_i = Int64[]            # first-stage solution
    y_i = Array{Int64,3}            # second-stage solution
    s_i = Array{Int64,2}            # second-stage slack variables
    it = 0              # iteration count
    
    scenario_based_model = build_scenario_based(inst, K)
    separation_model = build_separation_problem(K, inst)

    while (isempty(N) == false) && (runtime <= 240.0) # stop after 240 s
        it = it + 1
        # select unexplored node (TODO: which one to select?) and delete from set of unexplored nodes
        unique!(N)
        tau = popfirst!(N)
        # update model
        update_scenario_based!(scenario_based_model, tau)
        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        theta, x, y, s = solve_scenario_based_inplace(scenario_based_model, time_start)

        if theta < theta_i
            # update separation problem
            update_separation_problem!(separation_model, y, s)
            #(ζ, \xi, z)$ = Solve Separation Problem (8): max $ζ$ where $ζ$ is the amount of violation of the uncertain constraints and $\xi$ is the scenario that leads to violation
            zeta, xi = solve_separation_problem_inplace(separation_model, time_start)
            #println("separation problem solved, worst case scenario xi = $(xi)")

            if zeta <= 1e-6 # no violations
                # update incumbent
                theta_i = theta
                x_i = x
                y_i = y
                s_i = s
            else
                Knew = number_of_childnodes(tau)
                branch_partition!(N, tau, xi, Knew)
            end

        end
        runtime = (now()-time_start).value/1000
    end
    
    if theta_i == 1e10
        return "infeasible"
    else
        return x_i, y_i, s_i, theta_i, it, runtime
    end
end

"""
    create_scenario_based(inst, K)

Create the scenario-based K_adaptable problem for the uncertainty sets tau.
"""
function build_scenario_based(inst::AllocationInstance, K::Int64)

    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 2)
    J = size(loc_J, 2)
    W = inst.W
    D = inst.D

    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_bb_inplace); add_bridges = false)  # create model, reuse gurobi environment
    set_optimizer_attribute(rm, "OutputFlag", 0)    # no gurobi output
    set_string_names_on_creation(rm, false) # disable string names for performance improvement

    @expression(rm, c[i=1:I,j=1:J], norm(loc_I[:,i]-loc_J[:,j])); # transportation costs
    @expression(rm, slack_coeff, 10*max(c...))    # coefficient for slack variables in objective

    @variable(rm, 0 <= w[1:I] <= W, Int)            # first-stage decision
    @variable(rm, 0 <= q[1:I,1:J,1:K] <= W, Int)    # Second stage, wait-and-see decision how to distribute and slack
    @variable(rm, 0 <= s[1:J, 1:K] <= D, Int)       # One set of variables per cell
    
    @constraint(rm, sum(w[i] for i in 1:I) <= W)                          # supply limit
    @constraint(rm, [i=1:I,k=1:K], sum(q[i,j,k] for j in 1:J) <= w[i])    # service point limit

    @constraint(rm, demand_con[k=1:K, j=1:J], sum(q[i,j,k] for i in 1:I)+s[j,k] >= 0) # dummy demand constraint
        
    @variable(rm, 0 <= obj <= 10^10)                        # The objective function will be the maximum of the objective function across all the cells
    @variable(rm, 0<= z[1:K]<=10^10)                        # A variable to track each cells objective value
    @objective(rm, Min, obj)

    # Constrain objective function for this cell
    @constraint(rm, [k=1:K], z[k] >= slack_coeff*sum(s[j,k] for j in 1:J) + sum(c[i,j]*q[i,j,k] for i in 1:I, j in 1:J))
    @constraint(rm, [k=1:K], obj >= z[k])

    return rm
end

function update_scenario_based!(model::Model, tau::Vector{Vector{Vector{Int64}}}) 
    # Can't have constraints in place and just change rhs, 
    # because the number of constraints changes with tau and could be bigger or smaller than before.
    
    # delete old demand constraints
    delete.(model, model[:demand_con])
    unregister.(model, :demand_con)

    K = length(tau)
    J = size(model[:s],1)
    I = size(model[:w],1)
    # add new demand constraints
    @constraint(model, demand_con[k=1:K, xi=tau[k], j=1:J], sum(model[:q][i,j,k] for i in 1:I)+model[:s][j,k] >= xi[j])
    nothing #return model
end

function solve_scenario_based_inplace(model::Model, time::DateTime)
    set_remaining_time(model, time)
    # solve model
    optimize!(model)
    theta = objective_value(model)
    x = round.(Int, value.(model[:w]))
    y = round.(Int, value.(model[:q]))
    s = value.(model[:s])
    return theta, x, y, s
end


function build_separation_problem(K::Int64, inst::AllocationInstance)
    
    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 2)
    J = size(loc_J, 2)
    D = inst.D
    pc = inst.pc
    
    us = Model(() -> Gurobi.Optimizer(GRB_ENV_bb_inplace); add_bridges = false)
    set_optimizer_attribute(us, "OutputFlag", 0)
    set_string_names_on_creation(us, false) # disable string names for performance improvement

    @variable(us, zeta)     # amount of violation
    @variable(us, 0<= d[1:J] <=D)   # demand scenario // TODO: does it need to be Int?
    @variable(us, z[1:K,1:J], Bin)   # violation indicator
    @variable(us, y[1:I,1:J,1:K] >= 0) #y-variable to be fixed to the current solution
    #fix.(y, y_value; force = true)
    @variable(us, s[1:J,1:K] >= 0) #s-variable to be fixed to the current solution
    #fix.(s, s_value; force = true)

    # d must be in the uncertainty set
    @constraint(us, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))   # bound on aggregated demand
    for (j2,j1) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(us, d[j1]-d[j2] <= norm(loc_J[:,j1]-loc_J[:,j2],Inf))
    end

    @expression(us, M, 2*D+1)
    # only one (the worst) violated j for every k
    @constraint(us, [k=1:K], sum(z[k,j] for j in 1:J) == 1)
    # zeta measures the smallest violation among the k
    @constraint(us, [k=1:K, j=1:J], zeta + M*z[k,j] <= M + d[j] -( sum(y[i,j,k] for i in 1:I)+s[j,k]))

    @objective(us, Max, zeta)
    return us
end

function update_separation_problem!(sep_model::Model, y_value::Array{Int64, 3},s_value::Array{Int64,2})
    # fix variables y,s to current solution
    fix.(sep_model[:y], y_value; force = true)
    fix.(sep_model[:s], s_value; force = true)
    nothing #return sep_model
end

function solve_separation_problem_inplace(sep_model::Model, time::DateTime)
    set_remaining_time(sep_model, time)
    # optimize model
    optimize!(sep_model)
    return round.(value.(sep_model[:zeta]), digits = 4), ceil.(Int, value.(sep_model[:d]))#round.(value.(sep_model[:d]), digits = 2)
end