using JuMP, Gurobi
#const GRB_ENV_bb = Gurobi.Env()

"""
    solve_bb_general(K, inst)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. 
"""
function solve_bb_general(K, inst)
    # set of unexplored nodes, each node consists of K (disjoint) subsets of the uncertainty set
    # we start with K empty sets
    N = [[Iterators.repeated([],K)...]]
    # the incumbent
    theta_i = 10^10     # objective value
    x_i = []            # first-stage solution
    y_i = []            # second-stage solution
    s_i = []
    it = 0
    time_start = now()
    runtime = 0
    while (isempty(N) == false) && (runtime <= 240) # stop after 240 s
        it = it + 1
        #println("number of unexplored nodes: $(length(N))")
        # select unexplored node (TODO: which one to select?)
        # and delete from set of unexplored nodes
        
        tau = popfirst!(N)
        println("tau = $tau")
        #println("tau = $(tau)")
        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        @show theta, x, y, s = solve_scenario_based(tau, inst)
        #println("theta = $(theta), theta_i = $(theta_i)")
        if theta < theta_i

            #(ζ, \xi, z)$ = Solve Separation Problem (8): max $ζ$ where $ζ$ is the amount of violation of the uncertain constraints and $\xi$ is the scenario that leads to violation
            @show zeta, xi = solve_separation_problem_general(y, s, inst)
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
function solve_scenario_based(tau, inst)

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

    #rm = Model(() -> Gurobi.Optimizer(GRB_ENV_bb))
    rm = Model(Gurobi.Optimizer)
    set_optimizer_attribute(rm, "OutputFlag", 0)

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

function solve_separation_problem_general(y, s, inst)
    
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

"""
    solve_bb_tailored(K, inst)

Solve the K-adaptable problem with the Branch-and-Bound approach of Subramanyam et al. adapted to the supply pre-allocation model
"""
function solve_bb_tailored(K, inst)
    # set of unexplored nodes, each node consists of K (disjoint) subsets of the uncertainty set
    # we start with K empty sets
    N = [[Iterators.repeated([],K)...]]
    # the incumbent
    theta_i = 10^10     # objective value
    x_i = []            # first-stage solution
    y_i = []            # second-stage solution
    it = 0
    while (isempty(N) == false) && (it < 1000)
        it = it +1
        println("number of unexplored nodes: $(length(N))")
        println(N)
        # select unexplored node (TODO: which one to select?)
        # delete from set of unexplored nodes
        @show tau = popfirst!(N)

        # (θ, x, y) = Solve Scenario-based K-adapt Problem (6): min theta with uncsets tau 
        theta, x, y, s = solve_scenario_based(tau, inst)

        if theta < theta_i

            J = size(inst.loc_J, 1)  
            zeta = zeros(J)
            xi = zeros(J,J)
            for j_point in 1:J
                # Solve Separation Problem zeta is the amount of violation of the uncertain constraints and $\xi$ is the demand scenario that leads to violation
                zeta[j_point], xi[j_point,:] = solve_separation_problem_tailored(y, s, inst, j_point)
            end

            if all(<=(0), zeta) # no violations

                #$(θ^i, x^i, y^i) ← (θ, x, y)$
                theta_i = copy(theta)
                x_i = copy(x)
                y_i = copy(y)

            else
                zeta_bad = findall(z->z > 0, zeta) # the demand points that are violated
                N_bad = [copy(tau)] # the nodes/partitions to be added to the BB tree, initiate with tau (in the end we will have K^(J_bad) child nodes)

                for j in zeta_bad # for every violated demand site
                    iterationno = size(N_bad, 1) # for every partition in N_bad
                    for it in iterationno
                        tautau = popfirst!(N_bad) # pick out the partition
                        # Knew child nodes
                        # sort the partition by size of the subsets
                        sort!(tautau, by= x-> size(x), rev = true)
                        # define Knew the first empty subset
                        if isempty(tautau[end]) 
                            Knew = findfirst(x-> isempty(x), tautau)
                        else
                            Knew = K
                        end
                        for k in 1:Knew # add back a child node for every possibility to add the violation scenario of j

                            # each child node is the current uncset configuration...
                            tau_temp = copy(tautau)
                            # ...with the new scenario added to uncset number k
                            tau_temp[k] = union(tau_temp[k], [xi[j,:]])
                            push!(N_bad, copy(tau_temp))   
                        end
                    end
                end
                N = union(N, N_bad) # add all new partitions to N
            end

        end
    end
    
    if theta_i == 10^10

        return "infeasible"

    else

        return x_i, y_i, theta_i, it

    end
end

function solve_separation_problem_tailored(y, s, inst, j_point)
    
    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(inst.loc_I, 1)
    J = size(inst.loc_J, 1)  
    K = size(y, 3)
    D = inst.D  
    pc = inst.pc
        
    us = Model(Gurobi.Optimizer)

    @variable(us, zeta[1:K]) # max violation in case of plan k
    @variable(us, viola) # best-case violation
    @variable(us, 0<= d[1:J] <= D)   # demand scenario

    # d must be in the uncertainty set
    @constraint(us, sum(d[j] for j in 1:J) <= round(Int, pc*D*J))   # bound on aggregated demand
    for (j1,j2) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(us, d[j1]-d[j2] <= norm(loc_J[j1,:]-loc_J[j2,:],Inf))
    end

    @constraint(us, [k=1:K],  d[j_point] - (sum(y[i,j_point,k] for i in 1:I)+s[j_point,k])>= zeta[k] ) # demand violation at point j for plan k
    @constraint(us, [k=1:K], zeta[k] >= viola) # take the smallest violation (only one plan needs to be feasible)

    @objective(us, Max, viola)
    optimize!(us)

    return round.(value.(viola), digits = 4), round.(value.(d), digits = 2)
end
