using LinearAlgebra, Plots

include("helpers.jl")
include("solve_iter.jl")
include("solve_bb.jl")
include("solve_boxes.jl")
include("solve_bb_inplace.jl")
include("solve_boxes_inplace.jl")

number_of_instances = 20
times = zeros(number_of_instances, 5)
objectives = zeros(number_of_instances,5)
iterations = zeros(number_of_instances, 5)
no_sp = 1
no_dp = 5
k = 3
io = open("results/main.txt", "w")
write(io, "Solution times for $(no_sp) service points, $(no_dp) demand points, k=$(k)\norder: [BB, BB in place, box, box in place, iter]\n")
close(io)
rand_seeds = rand(500:999,number_of_instances)


for n in 1:number_of_instances
    
    se = rand_seeds[n]
    io = open("results/main.txt", "a")
    write(io, "seed: $se")
    close(io)
    #se=8

    inst_gen = generate_instance(no_sp,no_dp,se)
    #inst_gen = AllocationInstance([2 4; 8 3; 9 7], [1 9; 3 2; 4 1; 6 3; 7 9; 9 10], 15, 5, 0.5)
    #inst_gen = AllocationInstance([5 5], [5 2; 10 5], 5, 5, 0.5)

    # println("Locations:\nloc_I = $(inst_gen.loc_I)\nloc_J = $(inst_gen.loc_J)")

    println("BB started")
    x_general, y_general, s_general, theta_general, it_general, duration_general = solve_bb_general(k, inst_gen)

    println("BB in place started")
    x_general_ip, y_general_ip, s_general_ip, theta_general_ip, it_general_ip, duration_general_ip = solve_bb_inplace(k, inst_gen)

    println("box started")
    x_box, y_box, s_box, xi_box, theta_box, it_box, duration_box = solve_boxes(k, inst_gen)

    println("box in place started")
    x_box_ip, y_box_ip, s_box_ip, xi_box_ip, theta_box_ip, it_box_ip, duration_box_ip = solve_boxes_inplace(k, inst_gen)

    println("iter started")
    starting = now()
    theta_iter, x_iter, y_iter, k_iter, ktrue_iter = k_adapt_solution(2, inst_gen)
    ending = now()
    duration_iter = (ending-starting).value/1000 # convert to seconds

    times[n,:] = [duration_general duration_general_ip duration_box duration_box_ip duration_iter]

    objectives[n,:] = [theta_general theta_general_ip theta_box theta_box_ip theta_iter[end]]

    iterations[n,:] = [it_general it_general_ip it_box it_box_ip ktrue_iter[end]] # for iter, its the number of subsets not iterations

    #U = enum_uncset(demand_bound, loc_J_gen, cont_perc)
    #Obs = obs_obj(U, find_plan, q_v, c_gen)

    io = open("results/main.txt", "a")
    write(io, ",\t times: $(times[n,:]) s \t objectives: $(objectives[n,:]) \t iterations: $(iterations[n,:])\n")
    close(io)

end

# when was max runtime reached?
maxit_bb = findall(x-> x>240, times[:,1])
maxit_bb_ip = findall(x-> x>240, times[:,2])
maxit_box = findall(x-> x>240, times[:,3])
maxit_box_ip = findall(x-> x>240, times[:,4])


plot(1:number_of_instances, times, 
        label = ["BB" "BB_ip" "Box" "Box_ip" "iter"],
        ylabel="seconds", 
        xlabel="instance", 
        title="comp times for k=$k, $(number_of_instances) instances of size $(no_sp)x$(no_dp)")

# if length(maxit[:,1]) > 0
#     scatter!(maxit_bb, times[maxit_bb,1],label ="cutoff bb")
# end
# if length(maxit_box)> 0
#     scatter!(maxit_box, times[maxit_box,3],label ="cutoff box")
# end
savefig("results/times.png")

plot(1:number_of_instances, objectives, 
        label = ["BB" "BB_ip" "Box" "Box_ip" "iter"],
        ylabel="objective", 
        xlabel="instance", 
        title="objectives for k=$k, $(number_of_instances) instances of size $(no_sp)x$(no_dp)")

if length(maxit_bb) > 0
    scatter!(maxit_bb, objectives[maxit_bb,1], label ="cutoff bb")
end
if length(maxit_bb_ip) > 0
    scatter!(maxit_bb_ip, objectives[maxit_bb_ip,2], label ="cutoff bb_ip")
end
if length(maxit_box)> 0
    scatter!(maxit_box, objectives[maxit_box,3], label="cutoff box")
end
if length(maxit_box_ip)> 0
    scatter!(maxit_box_ip, objectives[maxit_box_ip,4], label="cutoff box_ip")
end
savefig("results/objectives.png")
#k_curve(o_v, p_v, n_val=no_sp, m_val=no_dp, rel_val=true) 


# k_curve_obs(p_v, Obs, n_val=no_sp, m_val=no_dp)

#obj_f(q, d, c) = 10*sum(max(d[j]-sum(q[i,j] for i in 1:no_sp),0) for j in 1:no_dp) + sum(c[i,j]*q[i,j] for i in 1:no_sp, j in 1:no_dp)
