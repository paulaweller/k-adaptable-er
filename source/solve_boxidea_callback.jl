using JuMP, Gurobi
const GRB_ENV_box = Gurobi.Env()

"""
    solve_boxes(K, loc_i, loc_J, W, D, pc)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_boxes(K::Int, inst::AllocationInstance)

    # if  

    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 1)
    J = size(loc_J, 1)
    c = reshape([norm(loc_I[i,:]-loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    D = inst.D
    pc = inst.pc
    T = length(tau)
    c = reshape([norm(loc_I[i,:]-loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    # coefficient for slack variables in objective
    slack_coeff = 10*max(c...)

    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_box))
    set_optimizer_attribute(rm, "OutputFlag", 0)

    # set solver time limit 
    set_time_limit_sec(rm, 240)

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

    # Callback checks for feasibility of partition
    function my_callback_function(cb_data)
        y_val = callback_value.(cb_data, q)
        s_val = callback_value.(cb_data, s)
        xi_val = callback_value.(cb_data, ξ)
        # solve separation problem
        zeta, u = solve_separation_problem_boxes(inst, K, xi_val)

        #push!(worst_case_scenarios, d_scenario)
        if zeta > 10^(-6) # a violation occurred
            con = @build_constraint(sum(q[i,j,p] for i in 1:I)+s[j,p] >= d_val)
            MOI.submit(rm, MOI.LazyConstraint(cb_data), con)
        end
    end
    MOI.set(rm, MOI.LazyConstraintCallback(), my_callback_function)


    @constraint(rm, [t=1:T], sum(v[t,k] for k=1:K) >= 1)        # every demand scenario must be covered by at least one plan
    @constraint(rm, [j=1:J,t=1:T,k=1:K], tau[t][j]*v[t,k] <= ξ[j,k])   # if plan k covers scenario tau[t], it must be componentwise larger

    # solve
    optimize!(rm)
    theta = getvalue(obj)
    x = round.(Int,getvalue.(w))
    y = round.(Int,getvalue.(q))
    s = round.(Int,getvalue.(s))
    xi = getvalue.(ξ)
    return x, y, s, xi, theta, runtime
end

function solve_separation_problem_boxes(inst::AllocationInstance, K::Int, ξ, time_start)

    loc_J = inst.loc_J
    J = size(loc_J, 1)
    D = inst.D
    pc = inst.pc
    us = Model(() -> Gurobi.Optimizer(GRB_ENV_box))
    set_optimizer_attribute(us, "OutputFlag", 0)

    @variable(us, zeta)     # amount of violation
    @variable(us, 0<= d[1:J] <=D)   # demand scenario // TODO: does it need to be Int?
    @variable(us, z[1:K,1:J], Bin)   # violation indicator

    # d must be in the uncertainty set
    @constraint(us, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))   # bound on aggregated demand
    for (j1,j2) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(us, d[j1]-d[j2] <= norm(loc_J[j1,:]-loc_J[j2,:],Inf))
    end

    M = 2*D+1
    @constraint(us, [k=1:K], sum(z[k,j] for j in 1:J) == 1)
    #@constraint(us, [k=1:K, j=1:J], zeta + ξ[j,k] <= d[j]*z[k,j])
    @constraint(us, [k=1:K, j=1:J], zeta + M*z[k,j] <= M + d[j] - ξ[j,k])

    @objective(us, Max, zeta)
    optimize!(us)

    return round.(value.(zeta), digits = 4), round.(value.(d), digits = 2)
end