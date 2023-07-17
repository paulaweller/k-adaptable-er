#-----------------------------------------------------------------------
# Apply the method from "Multistage Robust Mixed Integer Optimization
# with Adaptive Partitions" by Bertsimas and Dunning
#-----------------------------------------------------------------------

using JuMP, LinearAlgebra, Gurobi
const GRB_ENV_iter = Gurobi.Env()

import Random

"""
    TreeScenario
Stores the values of "active" uncertain parameters, as well as the
associated tree structure described in the paper.
"""
mutable struct TreeScenario
    ξ::Vector{Float64}
    parent
    children::Vector
end
is_leaf(t::TreeScenario) = isempty(t.children)

"""
    solve_partitioned_problem(loc_I, loc_J, W, D, pc, scenarios)
Solve the two-stage problem with one partition of the uncertainty
set for every leaf scenario in the `scenario_tree`. At optimality, grows the
tree given the new scenarios obtained.
"""
function solve_partitioned_problem(inst::AllocationInstance,
                                   scenario_tree::Vector{TreeScenario})

    I = size(inst.loc_I, 2)
    J = size(inst.loc_J, 2)

    # scenario_tree is a vector containing every scenario in the tree
    # We will have one cell for every leaf scenario in the tree
    leaf_scenarios = filter(is_leaf, scenario_tree)
    P = length(leaf_scenarios)  # Number of cells

    # Initialize the RO model
    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_iter); add_bridges = false)
    set_silent(rm)
    set_string_names_on_creation(rm, false) # disable string names for performance improvement
    set_optimizer_attribute(rm, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)

    @expression(rm, c[i=1:I,j=1:J], norm(inst.loc_I[:,i]-inst.loc_J[:,j])); # transportation costs
    @expression(rm, slack_coeff, 10.0*norm(c,Inf))                  # coefficient for slack variables in objective

    # Decision variables:
    # First stage, here-and-now decision where to store supplies
    @variable(rm, 0 <= w[1:I] <= inst.W, Int)
    # Second stage, wait-and-see decision how to distribute and slack
    # One set of variables per cell
    @variable(rm, 0 <= q[1:I,1:J,1:P] <= inst.W, Int)
    @variable(rm, 0 <= s[1:J, 1:P] <= inst.D, Int)
    # supply limit
    @constraint(rm, sum(w[i] for i in 1:I) <= inst.W)
    # The objective function will be the maximum of the objective function
    # across all the cells. Put a default upper bound, just so we don't
    # start off unbounded if we are using a cutting plane method.
    @variable(rm, 0 <= obj <= 10^10)
    # A variable to track each cells objective alue
    @variable(rm, 0<= z[1:P]<=10^10)
    # Minimize not only the maximum of the cells but also each cell objective
    @objective(rm, Min, obj+0.1*sum(z[i] for i in 1:P))

    # Constrain objective function for cells
    @constraint(rm, [p=1:P], z[p] >= slack_coeff*sum(s[j,p] for j in 1:J) + sum(c[i,j]*q[i,j,p] for i in 1:I, j in 1:J))
    @constraint(rm, [p=1:P], obj >= z[p])

    # service point limit
    @constraint(rm, [i=1:I, p=1:P], sum(q[i,j,p] for j in 1:J) <= w[i])

    # demand satisfaction as lazy callback
    # extract worst-case scenarios
    worst_case_scenarios = Vector{Float64}[]
    function my_callback_function(cb_data)
        y_val = callback_value.(cb_data, q)
        s_val = callback_value.(cb_data, s)
        for p in 1:P
            for j in 1:J
                d_val, d_scenario = solve_sep(p, j, inst.pc, inst.D, inst.loc_J, leaf_scenarios)
                push!(worst_case_scenarios, d_scenario)
                if sum(y_val[i,j,p] for i in 1:I)+s_val[j,p] < d_val
                    con = @build_constraint(sum(q[i,j,p] for i in 1:I)+s[j,p] >= d_val)
                    MOI.submit(rm, MOI.LazyConstraint(cb_data), con)
                end
            end
        end
    end
    MOI.set(rm, MOI.LazyConstraintCallback(), my_callback_function)

    # Solve
    optimize!(rm)                         
    # Extend the scenario tree
    for p in 1:P
        # if this cell's objective equals the worst cell's objective...
        if abs(value(z[p]) - value(obj)) < 1e-6
            for j in 1:J
                # Extract the active uncertain parameter values
                demand_scen = worst_case_scenarios[J*(p-1)+j]
                # Create a new child in the tree under this leaf
                demand_child = TreeScenario(demand_scen, leaf_scenarios[p], [])
                # Add to the tree
                push!(leaf_scenarios[p].children, demand_child)
                push!(scenario_tree, demand_child)
            end
        end
    end

    # Calculate actual number of plans
    q_original = round.(Int, value.(q))
    q_union = union(q_original[:,:,p] for p in 1:size(q_original, 3))
    n_plans = length(q_union)

    # Return the objective function value and the first-stage solution
    value(obj), value.(w), value.(q), P, n_plans
end


"""
    solve_sep(p, dn, pc, D, loc_J, scenario_tree)
Solve the separation problem for cell p and demand node dn. Returns worst-case d[dn] and d.
"""
function solve_sep(p::Int64, dn::Int64, pc::Float64, D::Float64, loc_J::Matrix{Float64}, scenario_tree)
    J = size(loc_J,2)
    leaf_scenarios = filter(is_leaf, scenario_tree)
    # Define the separation model
    sm = Model(() -> Gurobi.Optimizer(GRB_ENV_iter); add_bridges = false)
    set_silent(sm)
    set_string_names_on_creation(sm, false) # disable string names for performance improvement
    set_optimizer_attribute(sm, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)
    # variables
    @variable(sm, 0 <= d[1:J] <= D)
    # bound on aggregated demand
    @constraint(sm, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))

    # for each pair of demand points, add constraint that if locations are close, demand values must be close, too
    for j2 in 1:J
        for j1 in j2+1:J
            @constraint(sm, d[j1]-d[j2] <= norm(loc_J[:,j1]-loc_J[:,j2],Inf))
            @constraint(sm, d[j2]-d[j1] <= norm(loc_J[:,j1]-loc_J[:,j2],Inf))
        end
    end

    # We define multiple hyperplanes by walking up the scenario tree from this leaf.
    current_scenario = leaf_scenarios[p]
    parent_scenario = current_scenario.parent
    # We keep going until we hit the root of the tree, which is a scenario
    # that has no parent
    while parent_scenario !== nothing
        for sibling_scenario in parent_scenario.children
            if current_scenario == sibling_scenario
                continue  # Don't partition against ourself!
            end
            ξ_sub = sibling_scenario.ξ - current_scenario.ξ
            ξ_add = sibling_scenario.ξ + current_scenario.ξ
            @constraint(sm, dot(ξ_sub, d) <= dot(ξ_sub,ξ_add)/2)
        end
        # Move up the scenario tree
        current_scenario  = parent_scenario
        parent_scenario = current_scenario.parent
    end

    @objective(sm, Max, d[dn])
    optimize!(sm)
    worst_case_value = objective_value(sm)
    return worst_case_value, value.(d)
end

"""
    k_adapt_solution(it, inst)
Solve the problem for the parameters with it iterations.
"""
function k_adapt_solution(it::Int64, inst::AllocationInstance)
    # Start with no partitions (i.e., one scenario)
    scenario_tree = [ TreeScenario(zeros(4),nothing,[]) ]
    # store these values for every iteration:
    obj_val = Float64[]    # objective value
    w_val = Vector{Float64}[]    # storage of supplies
    q_val = Vector{Array{Float64,3}}[]      # transport plans
    p_val = Int64[]      # number pf cells
    p_true = Int64[]     # number of plans

    # q_it = 0
    for i in 1:it
        #println("Iteration $i started.")
        obj_it, w_it, q_it, p_it, p_true_it = solve_partitioned_problem(inst, scenario_tree)
        push!(obj_val, obj_it)
        push!(w_val, w_it)
        push!(q_val, [q_it])
        push!(p_val, p_it)
        push!(p_true, p_true_it)
    end    

    #= println("k-adaptable solution")
    println("objectives: $obj_val")
    println("supply status = $w_val")
    println("transportation: $(q_val)")
    println("number of cells: $p_val")
    println("actual number of cells: $p_val") =#
    return obj_val, w_val, q_val, p_val, p_true
end



