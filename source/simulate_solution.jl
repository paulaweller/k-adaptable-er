include("helpers_data.jl")

using JuMP, Gurobi, BilevelJuMP
const GRB_ENV_bb = Gurobi.Env()

function observable_worst_case_objectives(no, mo, pco, ko)
    pb, bb, box = extract_solutions(no, mo, pco, ko)

    obj_pb = []
    obj_bb = []
    obj_box = []
    objo_pb = []
    objo_box = []

    for (key, values) in pb
        # Extract x and y coordinates from each 2D vector
        
        val_bb = bb[key]
        val_box = box[key]

        # key is [n,m,pc,k,l]   
        n = key[1]
        m = key[2]
        pc = key[3]
        k = key[4]
        l = key[5]   
        it = round(Int, 1 + (k-1)/(m-1))  
        if val_bb == "n"
            
            y_val_pb = arrayfromstr(values, n,m,it)
            y_val_box = arrayfromstr(val_box, n,m,k)
            # read problem instance
            inst = read_instance_for_param(n,m,pc,l)

            theta_pb = solve_worst_case_objective(y_val_pb, inst)
            theta_box = solve_worst_case_objective(y_val_box, inst)

            push!(objo_pb, theta_pb)
            push!(objo_box, theta_box)
        elseif val_bb =='n'  

            y_val_pb = arrayfromstr(values, n,m,it)
            y_val_box = arrayfromstr(val_box, n,m,k)
            # read problem instance
            inst = read_instance_for_param(n,m,pc,l)

            theta_pb = solve_worst_case_objective(y_val_pb, inst)
            theta_box = solve_worst_case_objective(y_val_box, inst)

            push!(objo_pb, theta_pb)
            push!(objo_box, theta_box)
        else  

            y_val_pb = arrayfromstr(values, n,m,it)
            y_val_bb = arrayfromstr(val_bb, n,m,k)
            y_val_box = arrayfromstr(val_box, n,m,k)
            # read problem instance
            inst = read_instance_for_param(n,m,pc,l)

            theta_pb = solve_worst_case_objective(y_val_pb, inst)
            theta_bb = solve_worst_case_objective(y_val_bb, inst)
            theta_box = solve_worst_case_objective(y_val_box, inst)

            push!(obj_pb, theta_pb)
            push!(obj_bb, theta_bb)
            push!(obj_box, theta_box)
        end
        
    end
    return obj_pb, obj_bb, obj_box, objo_pb, objo_box 
end

function arrayfromstr(str, n, m, k)
    final = zeros(convert(Int,n),convert(Int,m),convert(Int,k))
    # split by k
    str = split(chop(str, head=1), ";;;")
    str = strip.(str, ['['])
    str = strip.(str, [']'])
    k_it = 0
    for str_k in str
        if length(str_k) > 0
            k_it = k_it+1
            # remove spaces in beginning and end
            str_k = strip(str_k, [' '])
            # split by n (service points)
            str_k = split(str_k, ";")
            n_it = 0
            for str_k_n in str_k
                n_it = n_it+1
                str_k_n = strip(str_k_n, [' '])
                str_k_n = split(str_k_n, " ")
                if length(str_k_n) == m
                    str_k_n = parse.(Float64,str_k_n)
                    final[n_it,:,k_it] = str_k_n
                #elseif str_k_n == [""]
                else error("Expected m = $m, but received m = $(length(str_k_n))")
                end
            end
            (n_it != n) ? error("Expected n = $n, but received n = $(n_it)") : nothing
        end
    end
    (k_it != k) ? error("Expected k = $k, but received k = $(k_it) for n = $n, m=$m") : nothing
    return final
end


"""
    extract_solutions(no, mo, pco, ko)

Returns a dictionary of the y-values of all solutions available for no service points, mo demand points, pco percentages, ko k's.
"""
function extract_solutions(no, mo, pco, ko)
    y_pb = Dict()
    y_bb = Dict()
    y_box = Dict()
    for n in no
        for m in mo
            for pc in pco
                for k in ko
                    # Specify the directory path
                    directory_path = "source/results/data_batch_$(n)_$(m)_$(pc)/individual"

                    # Specify the common beginning of the filename
                    common_prefix = "sol_data_batch_$(n)_$(m)_$(pc)_k$(k)_"

                    # Get a list of files in the directory that match the pattern
                    files = filter(x -> occursin(common_prefix, x), readdir(directory_path))
                    for file in files
                        l = file[(end-5):(end-4)]
                        l = strip(l, ['l'])
                        l = parse(Int, l)
                        dict = read_solutions_from_file("source/results/data_batch_$(n)_$(m)_$(pc)/individual/$(file)")

                        dict_pb = Dict([n,m,pc,k,l] => dict["y_pb"])
                        dict_bb = Dict([n,m,pc,k,l] => dict["y_bb"])
                        dict_box = Dict([n,m,pc,k,l] => dict["y_box"])
                        merge!(y_pb, dict_pb)
                        merge!(y_bb, dict_bb)
                        merge!(y_box, dict_box)
                    end
                end
            end
        end
    end
    return y_pb, y_bb, y_box
end

function read_instance_for_param(n,m,pc,l)
    n = convert(Int,n)
    m=convert(Int, m)
    file = "source/data/data_batch_$(n)_$(m)_$(pc).txt"
    inst = read_one_instance_from_file(convert(Int, l), file)
    return inst
end

"""
    solve_worst_case_objective(y, inst)

Solve the scenario-based K_adaptable problem for the uncertainty sets tau.
"""
function solve_worst_case_objective(y, inst::AllocationInstance)

    loc_I = inst.loc_I
    loc_J = inst.loc_J
    I = size(loc_I, 2)
    J = size(loc_J, 2)
    W = inst.W
    D = inst.D
    pc = inst.pc
    D_max = max(D...)
    K = size(y, 3)

    m = BilevelModel(() -> Gurobi.Optimizer(GRB_ENV_bb))
    set_optimizer_attribute(m, "OutputFlag", 0)

    ########### upper model

    @variable(Upper(m), 0 <= d[1:J])
    # d must be in the uncertainty set
    @constraint(Upper(m), [j=1:J], d[j] <= D[j])
    @constraint(Upper(m), sum(d[j] for j in 1:J) <= floor(pc*sum(D)))   # bound on aggregated demand
    for (j2,j1) in Iterators.product(1:J,1:J)   # clustering of demand
        @constraint(Upper(m), d[j1]-d[j2] <= norm(loc_J[:,j1]-loc_J[:,j2],Inf))
    end

    # available second-stage policies
    @variable(Upper(m), 0 <= q[1:I,1:J,1:K] <= W, Int)
    fix.(m[:q], y; force = true)



    ########## lower model

    @expression(Lower(m), c[i=1:I,j=1:J], norm(loc_I[:,i]-loc_J[:,j])); # transportation costs
    @expression(Lower(m), slack_coeff, 1000.0*norm(c,Inf))              # coefficient for slack variables in objective


    @variable(Lower(m), 0 <= obj <= 1e10)                        # The objective function will be the minimum of the objective function across all the policies
    @variable(Lower(m), 0 <= z[1:K] <= 1e10)                     # A variable to track each policies objective value
    @variable(Lower(m), 0 <= s[1:J, 1:K]<= D_max)                # Unsatisfied demand

    # Constrain objective function for each policy
    @constraint(Lower(m),[k=1:K], z[k] == slack_coeff*sum(s[j,k] for j in 1:J) + sum(c[i,j]*q[i,j,k] for i in 1:I, j in 1:J))
    # Final objective is the most beneficial policy
    @constraint(Lower(m), [k=1:K], obj <= z[k])

    # Demand must be satisfied
    @constraint(Lower(m), [j=1:J, k=1:K], s[j,k] >= d[j] - sum(q[i,j,k] for i in 1:I))

    # minimize cost for all policies first, then maximize obj to be the smallest one
    @objective(Lower(m), Min, sum(z[k] for k in 1:K) - 0.01*obj)
    @objective(Upper(m), Max, obj)

    # set_remaining_time(m, time_start, time_limit)
    # solve
    optimize!(m)
    # if result_count(rm) == 0
    #     return 1e10, zeros(Float64, I), zeros(Float64, I, J, K), zeros(Float64, J, K)
    # end
    theta = objective_value(m)

    return theta
end