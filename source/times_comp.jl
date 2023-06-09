using LinearAlgebra, Plots

include("helpers.jl")
include("solve_iter.jl")
include("solve_bb.jl")
include("solve_boxidea.jl")

number_of_instances = 25
times = zeros(number_of_instances, 3)
objectives = zeros(number_of_instances,3)
iterations = zeros(number_of_instances, 3)
no_sp = 2
no_dp = 3
k = 4
io = open("results/main.txt", "w")
write(io, "Solution times for $(no_sp) service points, $(no_dp) demand points, k=$(k)\norder: [BB, iter, box]\n")
close(io)
rand_seeds = rand(500:999,25)


for n in 1:number_of_instances
    
    se = rand_seeds[n]
    io = open("results/main.txt", "a")
    write(io, "seed: $se")
    close(io)
    #se=8

    inst_gen = generate_instance(no_sp,no_dp,se)
    #inst_gen = AllocationInstance([2 4; 8 3; 9 7], [1 9; 3 2; 4 1; 6 3; 7 9; 9 10], 15, 5, 0.5)
    #inst_gen = AllocationInstance([5 5], [5 2; 10 5], 5, 5, 0.5)

    I = size(inst_gen.loc_I, 1)
    J = size(inst_gen.loc_J, 1)
    c = reshape([norm(inst_gen.loc_I[i,:]-inst_gen.loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
    # println("Locations:\nloc_I = $(inst_gen.loc_I)\nloc_J = $(inst_gen.loc_J)")

    x_general, y_general, s_general, theta_general, it_general, duration_general = solve_bb_general(k, inst_gen)

    starting = now()
    theta_iter, x_iter, y_iter, k_iter, ktrue_iter = k_adapt_solution(2, inst_gen)
    ending = now()
    duration_iter = (ending-starting).value/1000 # convert to seconds

    x_box, y_box, s_box, xi_box, theta_box, it_box, duration_box = solve_boxes(k, inst_gen.loc_I, inst_gen.loc_J, inst_gen.W, inst_gen.D, inst_gen.pc)

    

    times[n, 1] = duration_general
    times[n, 2] = duration_iter
    times[n, 3] = duration_box

    objectives[n, 1] = theta_general
    objectives[n, 2] = theta_iter[end]
    objectives[n, 3] = theta_box

    iterations[n, 1] = it_general
    iterations[n, 2] = ktrue_iter[end]   # number of partitions
    iterations[n, 3] = it_box

    #U = enum_uncset(demand_bound, loc_J_gen, cont_perc)
    #Obs = obs_obj(U, find_plan, q_v, c_gen)

    io = open("results/main.txt", "a")
    write(io, ",\t times: $(times[n,:]) s \t objectives: $(objectives[n,:]) \t iterations: $(iterations[n,:])\n")
    close(io)

end

# when was max runtime reached?
maxit_bb = findall(x-> x>240, times[:,1])
maxit_box = findall(x-> x>240, times[:,3])

plot(1:number_of_instances, times[:,1], ylabel="seconds", xlabel="instance", title="comp times for k=3, $(number_of_instances) instances of size $(no_sp)x$(no_dp)", label="BB")
plot!(1:number_of_instances, times[:,2], label="iter")
plot!(1:number_of_instances, times[:,3], label="Box")
if length(maxit_bb) > 0
    scatter!(maxit_bb, times[maxit_bb,1],label ="cutoff bb")
end
if length(maxit_box)> 0
    scatter!(maxit_box, times[maxit_box,3],label ="cutoff box")
end
savefig("results/times.png")

plot(1:number_of_instances, objectives[:,1], ylabel="objective", xlabel="instance", title="objectives for k=3, $(number_of_instances) instances of size $(no_sp)x$(no_dp)", label="BB")
plot!(1:number_of_instances, objectives[:,2], label="iter")
plot!(1:number_of_instances, objectives[:,3], label="Box")
if length(maxit_bb) > 0
    scatter!(maxit_bb, objectives[maxit_bb,1], label ="cutoff bb")
end
if length(maxit_box)> 0
    scatter!(maxit_box, objectives[maxit_box,3], label="cutoff box")
end
savefig("results/objectives.png")
#k_curve(o_v, p_v, n_val=no_sp, m_val=no_dp, rel_val=true) 


# k_curve_obs(p_v, Obs, n_val=no_sp, m_val=no_dp)

#obj_f(q, d, c) = 10*sum(max(d[j]-sum(q[i,j] for i in 1:no_sp),0) for j in 1:no_dp) + sum(c[i,j]*q[i,j] for i in 1:no_sp, j in 1:no_dp)
