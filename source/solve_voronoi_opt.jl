using JuMP, Ipopt, LinearAlgebra, Dualization, Gurobi
include("helpers.jl")


"""
    solve_voronoi_kkt(K, inst)

    Solve the K-adaptable pre-allocation problem for instance "inst" with Voronoi-style partitioning.
"""

function solve_voronoi_opt(K, inst)

    I_size = size(inst.loc_I, 1)
    J_size = size(inst.loc_J, 1)
    c = reshape([norm(inst.loc_I[i,:]-inst.loc_J[j,:]) for j in 1:J_size for i in 1:I_size],I_size,J_size)
    

    part_model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(part_model, "NonConvex", 2)
    # set_optimizer_attribute(part_model, "FeasibilityTol", 1e-4)
    set_optimizer_attribute(part_model, "TimeLimit", 120)
    set_optimizer_attribute(part_model, "Presolve", 0)
    # variables ---------------------------------

    # supply
    @variable(part_model, 0<=x[1:I_size]<=inst.W)
    # allocation of supplies
    @variable(part_model, 0<=y[1:I_size,1:J_size,1:K]<=inst.W)
    # unsatisfied demand
    @variable(part_model, 0<=s[1:J_size,1:K]<=inst.D)
    # objectives of cells
    @variable(part_model, z[1:K]>= 0)
    # overall objective
    @variable(part_model, z_obj>= 0)
    # voronoi points
    @variable(part_model, 0<=v[1:J_size,1:K])

    # voronoi points can't be equal
    @constraint(part_model, [k=1:K, l=vcat([1:k-1]...,[k+1:K]...)], sum( (v[i,k]-v[i,l])^2 for i in 1:J_size)>=0.1)

    # uncertainty variables
    # worst-case demand scenario, one per partition k and row j
    @variable(part_model, 0<= ξ[1:J_size,1:K]<=inst.D)


    # upper-level constraints ----------------------

    # z_obj is worst of the K cell objectives
    @constraint(part_model, [k=1:K], z[k] <= z_obj)
    # define cell objectives
    @constraint(part_model, [k=1:K], z[k] >= 100*sum(s[j,k] for j=1:J_size)+sum(c[i,j]*y[i,j,k] for i=1:I_size for j=1:J_size))
    # allocation must comply with supply storage
    @constraint(part_model, [i=1:I_size,k=1:K], sum(y[i,j,k] for j=1:J_size) <=x[i])
    # bounded supply
    @constraint(part_model, sum(x[i] for i=1:I_size)<= inst.W)
    # demand must be satisfied
    @constraint(part_model, [j=1:J_size, k=1:K], sum(y[i,j,k] for i=1:I_size)+s[j,k] >= ξ[j,k])

    # voronoi points can't coincide
    # TODO @constraint(part_model, [k=1:K, l=(k+1):K], ones(J_size)*(v[:,k] .- v[:,l]))

    # lower-level model --------------------------

    # dual variable for distance constraints for every j, k, j1, j2
    @variable(part_model, α[1:J_size,1:K,1:J_size,1:J_size]>=0)

    # dual variable for demand sum constraint
    @variable(part_model, β[1:J_size,1:K]>=0)

    # dual variable for voronoi constraints for every j,k,l
    @variable(part_model, γ[1:J_size, 1:K, 1:K]>=0)

    # dual objective equals ξ[j,k]
    dist = reshape([norm(inst.loc_J[j1,:]-inst.loc_J[j2,:],Inf) for j1 in 1:J_size for j2 in 1:J_size],J_size,J_size)
    @variable(part_model, v_dot[1:K,1:K]) #k,l
    @constraint(part_model, [k=1:K, l=vcat([1:k-1]...,[k+1:K]...)], v_dot[k,l] == sum((v[j,l]-v[j,k])*(v[j,l]+v[j,k]) for j in 1:J_size))

    @constraint(part_model, [j=1:J_size,k=1:K], sum(dist[j1,j2]*α[j,k,j1,j2] for j1 in 1:J_size for j2 in vcat([1:j1-1]...,[j1+1:J_size]...))
                                            + β[j,k]*inst.pc*inst.D*J_size + sum(0.5*v_dot[k,l]*γ[j,k,l] for l in vcat([1:k-1]...,[k+1:K]...)) <= ξ[j,k])

    @constraint(part_model, [j=1:J_size, k=1:K, j1=vcat([1:j-1]...,[j+1:J_size]...)], sum(α[j,k,j1,j2]-α[j,k,j2,j1] for j2 in vcat([1:j1-1]...,[j1+1:J_size]...)) 
                                                                                +β[j,k] +sum((v[j1,l]-v[j1,k])*γ[j,k,l] for l in vcat([1:k-1]...,[k+1:K]...))>=0)

    @constraint(part_model, [j=1:J_size, k=1:K], sum(α[j,k,j,j2]-α[j,k,j2,j] for j2 in vcat([1:j-1]...,[j+1:J_size]...)) 
    +β[j,k] +sum((v[j,l]-v[j,k])*γ[j,k,l] for l in vcat([1:k-1]...,[k+1:K]...))>=1)
    
    @objective(part_model, Min, z_obj)

    open("results/model.txt","a") do io
        println(io,part_model)
     end
    optimize!(part_model)

    return value(z_obj), value.(v), value.(x), value.(y), value.(s), value.(ξ)
end

