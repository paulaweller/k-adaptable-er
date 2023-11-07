using JuMP, Gurobi
const GRB_ENV_bb = Gurobi.Env()

"""
    solve_bb_general(K, inst)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_bb(K::Int64, inst::AllocationInstance; time_limit::Float64 = 240.0)
    time_start = now()
    runtime = 0.0
    I = size(inst.loc_I, 2)
    J = size(inst.loc_J, 2)
    # set of unexplored nodes, each node consists of K (disjoint) subsets of the uncertainty set
    # we start with K empty sets
    N = Vector{Vector{Vector{Float64}}}[]
    push!(N, Vector{Vector{Float64}}[Iterators.repeated(Vector{Float64}[],K)...])
    # the incumbent
    theta_i = 1e10     # objective value
    x_i = Float64[]            # first-stage solution
    y_i = zeros(Float64,I,J,K)            # second-stage solution
    s_i = zeros(Float64, J,K)           # second-stage slack variables
    best_partition = Vector{Vector{Float64}}[Iterators.repeated(Vector{Float64}[],K)...]
    obj_evolution = Vector{Float64}[]
    it = 0              # iteration count
    
    while (isempty(N) == false) && (runtime <= time_limit) # stop after 120 s
        it = it + 1
        #println("number of unexplored nodes: $(length(N))")
        # select unexplored node (TODO: which one to select?)
        # and delete from set of unexplored nodes
        
        tau = popfirst!(N)

        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        theta, x, y, s = solve_scenario_based(tau, inst, time_start, time_limit)
        #println("theta = $(theta), theta_i = $(theta_i)")
        if theta < theta_i

            #(ζ, \xi, z)$ = Solve Separation Problem (8): max $ζ$ where $ζ$ is the amount of violation of the uncertain constraints and $\xi$ is the scenario that leads to violation
            zeta, xi = solve_separation_problem_general(y, s, inst, time_start, time_limit)
            #println("separation problem solved, worst case scenario xi = $(xi)")

            if zeta <= 1e-3 # no violations

                #$(θ^i, x^i, y^i) ← (θ, x, y)$
                theta_i = theta
                x_i = x
                y_i = y
                s_i = s
                best_partition = tau
                #println("incumbent found at time ", 
                timenow = (now()-time_start).value/1000
                push!(obj_evolution, [copy(timenow), copy(theta)])

            else
                Knew = number_of_childnodes(tau)
                branch_partition!(N, tau, ceil.(round.(xi,digits = 4)), Knew)
                unique!(N)
            end

        end
        runtime = (now()-time_start).value/1000
    end
    
    if theta_i == 1e10

        return "infeasible"

    else

        return x_i, y_i, s_i, best_partition, theta_i, it, runtime, obj_evolution

    end
end

"""
    solve_scenario_based(tau, inst)

Solve the scenario-based K_adaptable problem for the uncertainty sets tau.
"""
function solve_scenario_based(tau::Vector{Vector{Vector{Float64}}}, inst::AllocationInstance, time_start::DateTime, time_limit::Float64)

    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 2)
    J = size(loc_J, 2)
    K = length(tau)
    W = inst.W
    D = inst.D
    D_max = max(D...)

    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_bb); add_bridges = false)
    set_optimizer_attribute(rm, "OutputFlag", 0)
    set_string_names_on_creation(rm, false) # disable string names for performance improvement
    set_optimizer_attribute(rm, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)
    set_optimizer_attribute(rm, "Threads", 8)

    @expression(rm, c[i=1:I,j=1:J], norm(loc_I[:,i]-loc_J[:,j])); # transportation costs
    @expression(rm, slack_coeff, 1000.0*norm(c,Inf))                  # coefficient for slack variables in objective

    @variable(rm, 0 <= w[1:I] <= W, Int)            # first-stage decision
    @variable(rm, 0 <= q[1:I,1:J,1:K] <= W, Int)    # Second stage, wait-and-see decision how to distribute and slack
    @variable(rm, 0 <= s[1:J, 1:K]<= D_max, Int)       # One set of variables per cell
    
    
    @constraint(rm, sum(w[i] for i in 1:I) <= W)                    # supply limit
    @constraint(rm, [i=1:I,k=1:K], sum(q[i,j,k] for j in 1:J) <= w[i])    # service point limit
        
    @variable(rm, 0 <= obj <= 1e10)                        # The objective function will be the maximum of the objective function across all the cells
    @variable(rm, 0<= z[1:K]<=1e10)                        # A variable to track each cells objective value
    @objective(rm, Min, obj)

    for k in 1:K
        
        # Constrain objective function for this cell
        @constraint(rm, z[k] >= slack_coeff*sum(s[j,k] for j in 1:J) + sum(c[i,j]*q[i,j,k] for i in 1:I, j in 1:J))
        @constraint(rm, obj >= z[k])

        # Demand must be satisfied
        for xi in tau[k]
            #println(xi)
            @constraint(rm, [j=1:J], sum(q[i,j,k] for i in 1:I)+s[j,k] >= xi[j])
        end
 
    end
    set_remaining_time(rm, time_start, time_limit)
    # solve
    optimize!(rm)
    if result_count(rm) == 0
        return 1e10, zeros(Float64, I), zeros(Float64, I, J, K), zeros(Float64, J, K)
    end
    theta::Float64 = objective_value(rm)
    x = value.(w)
    y = value.(q)
    s_val = value.(s)

    return theta, x, y, s_val
end

function solve_separation_problem_general(y::Array{Float64,3}, s::Array{Float64,2}, inst::AllocationInstance, time_start::DateTime, time_limit::Float64)
    
    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 2)
    J = size(loc_J, 2)
    K = size(y, 3)
    D = inst.D
    D_max = max(D...)
    pc = inst.pc
    
    us = Model(() -> Gurobi.Optimizer(GRB_ENV_bb); add_bridges = false)
    set_optimizer_attribute(us, "OutputFlag", 0)
    set_string_names_on_creation(us, false) # disable string names for performance improvement
    set_optimizer_attribute(us, "MIPGap", 1e-3) # set gap to 0.1% (default is 1e-4)
    set_optimizer_attribute(us, "Threads", 8)

    @variable(us, zeta)     # amount of violation
    @variable(us, 0 <= d[1:J])   # demand scenario // TODO: does it need to be Int?
    @constraint(us, [j=1:J], d[j] <= D[j])
    @variable(us, z[1:K,1:J], Bin)   # violation indicator

    # d must be in the uncertainty set
    @constraint(us, sum(d[j] for j in 1:J) <= floor(pc*sum(D)))   # bound on aggregated demand
    for (j2,j1) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(us, d[j1]-d[j2] <= norm(loc_J[:,j1]-loc_J[:,j2],Inf))
    end

    @expression(us, M, 2*D_max+1)
    @constraint(us, [k=1:K], sum(z[k,j] for j in 1:J) == 1)
    #@constraint(us, [k=1:K], zeta + M*z[k,0] <= M + slack_coeff*sum(s[j,k] for j in 1:J) + sum(c[i,j]*y[i,j,k] for i in 1:I, j in 1:J) - theta) #don't need this because objective coeffs are certain
    @constraint(us, [k=1:K, j=1:J], zeta + M*z[k,j] <= M + d[j] -( sum(y[i,j,k] for i in 1:I)+s[j,k]))

    @objective(us, Max, zeta)

    set_remaining_time(us, time_start, time_limit)
    optimize!(us)

    if result_count(us) == 0
        return 1, zeros(Float64, J)
    end

    return value.(zeta), value.(d)
end