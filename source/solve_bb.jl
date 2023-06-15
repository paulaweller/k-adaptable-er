using JuMP, Gurobi
const GRB_ENV_bb = Gurobi.Env()

"""
    solve_bb_general(K, inst)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_bb_general(K, inst)
    time_start = now()
    runtime = 0
    # set of unexplored nodes, each node consists of K (disjoint) subsets of the uncertainty set
    # we start with K empty sets
    N = [[Iterators.repeated([],K)...]]
    # the incumbent
    theta_i = 10^10     # objective value
    x_i = []            # first-stage solution
    y_i = []            # second-stage solution
    s_i = []
    it = 0
    
    while (isempty(N) == false) && (runtime <= 240) # stop after 240 s
        it = it + 1
        #println("number of unexplored nodes: $(length(N))")
        # select unexplored node (TODO: which one to select?)
        # and delete from set of unexplored nodes
        
        tau = popfirst!(N)

        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        theta, x, y, s = solve_scenario_based(tau, inst, time_start)
        #println("theta = $(theta), theta_i = $(theta_i)")
        if theta < theta_i

            #(ζ, \xi, z)$ = Solve Separation Problem (8): max $ζ$ where $ζ$ is the amount of violation of the uncertain constraints and $\xi$ is the scenario that leads to violation
            zeta, xi = solve_separation_problem_general(y, s, inst, time_start)
            #println("separation problem solved, worst case scenario xi = $(xi)")

            if zeta <= 10^(-6) # no violations

                #$(θ^i, x^i, y^i) ← (θ, x, y)$
                theta_i = copy(theta)
                x_i = copy(x)
                y_i = copy(y)
                s_i = copy(s)

            else
                # Knew child nodes
                # sort the partition by size of the subsets
                sort!(tau, by= x-> size(x), rev = true)
                #define Knew the first empty subset
                if isempty(tau[end]) 
                    Knew = findfirst(x-> isempty(x), tau)
                else
                    Knew = K
                end
                for k in 1:Knew
                    # each child node is the current uncset configuration...
                    tau_temp = copy(tau)
                    # ...with the new scenario added to uncset number k
                    tau_temp[k] = union(tau_temp[k], [xi])
                    #push!(tau_temp[k], xi)
                    N = union(N, [copy(tau_temp)])
                    
                    #push!(N, copy(tau_temp))   
                end
                #println("N updated, N = $N")

            end

        end
        runtime = (now()-time_start).value/1000
    end
    
    if theta_i == 10^10

        return "infeasible"

    else

        return x_i, y_i, s_i, theta_i, it, runtime

    end
end

"""
    solve_scenario_based(tau, inst)

Solve the scenario-based K_adaptable problem for the uncertainty sets tau.
"""
function solve_scenario_based(tau, inst, time_start)

    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 1)
    J = size(loc_J, 1)
    K = length(tau)
    c = reshape([norm(loc_I[i,:]-loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    W = inst.W
    D = inst.D
    # coefficient for slack variables in objective
    slack_coeff = 10*max(c...)

    rm = Model(() -> Gurobi.Optimizer(GRB_ENV_bb))
    set_optimizer_attribute(rm, "OutputFlag", 0)
    # calculate remaining time before cutoff
    time_remaining = 240 + (time_start - now()).value/1000
    # set solver time limit accordingly
    set_time_limit_sec(rm, time_remaining)

    @variable(rm, 0 <= w[1:I] <= W, Int)            # first-stage decision
    @variable(rm, 0 <= q[1:I,1:J,1:K] <= W, Int)    # Second stage, wait-and-see decision how to distribute and slack
    @variable(rm, 0 <= s[1:J, 1:K] <= D, Int)       # One set of variables per cell
    
    
    @constraint(rm, sum(w[i] for i in 1:I) <= W)                    # supply limit
    @constraint(rm, [i=1:I,k=1:K], sum(q[i,j,k] for j in 1:J) <= w[i])    # service point limit
        
    @variable(rm, 0 <= obj <= 10^10)                        # The objective function will be the maximum of the objective function across all the cells
    @variable(rm, 0<= z[1:K]<=10^10)                        # A variable to track each cells objective value
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

    # solve
    optimize!(rm)
    theta = objective_value(rm)
    x = round.(Int, value.(w))
    y = round.(Int, value.(q))
    s = value.(s)
    return theta, x, y, s
end

function solve_separation_problem_general(y, s, inst, time_start)
    
    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 1)
    J = size(loc_J, 1)
    K = size(y, 3)
    c = reshape([norm(loc_I[i,:]-loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    D = inst.D
    pc = inst.pc
    
    us = Model(() -> Gurobi.Optimizer(GRB_ENV_bb))
    set_optimizer_attribute(us, "OutputFlag", 0)
    # calculate remaining time before cutoff
    time_remaining = 240 + (time_start - now()).value/1000
    # set solver time limit accordingly
    set_time_limit_sec(us, time_remaining)

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
    #@constraint(us, [k=1:K], zeta + M*z[k,0] <= M + slack_coeff*sum(s[j,k] for j in 1:J) + sum(c[i,j]*y[i,j,k] for i in 1:I, j in 1:J) - theta) #don't need this because objective coeffs are certain
    @constraint(us, [k=1:K, j=1:J], zeta + M*z[k,j] <= M + d[j] -( sum(y[i,j,k] for i in 1:I)+s[j,k]))

    @objective(us, Max, zeta)
    optimize!(us)

    return round.(value.(zeta), digits = 4), round.(value.(d), digits = 2)
end