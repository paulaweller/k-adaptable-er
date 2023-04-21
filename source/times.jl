using LinearAlgebra, Plots, PrettyTables, Statistics

include("solve_static.jl")
include("solve_comp.jl")
include("solve_iter.jl")
include("helpers.jl")
times = zeros(25, 10)
io = open("results/main.txt", "w")
write(io, "K-adapt times for 5 service points, n demand points, 3 iterations\n")
close(io)
rand_seeds = rand(500:999, 10,25)

for count in 1:10
    io = open("results/main.txt", "a")
    write(io, "Iteration $count:\n")
    close(io)
    for n in 1:25
        no_sp = 5
        no_dp = n
        se = rand_seeds[count,n]
        io = open("results/main.txt", "a")
        write(io, "service points: $no_sp,\t demand points: $no_dp,\t seed: $se")
        close(io)
        #se=8

        loc_I_gen, loc_J_gen, demand_bound, cont_perc, agg_supply_bound = generate_instance(no_sp,no_dp, se, plot_loc=false, loc_max=20)
        #write(io, "loc_I_gen = $loc_I_gen, loc_J_gen = $loc_J_gen, demand_bound = $demand_bound, cont_perc = $cont_perc \n")

        #doc example
        # loc_I_gen = [2 4; 8 3; 9 7]
        # loc_J_gen = [1 9; 3 2; 4 1; 6 3; 7 9; 9 10]
        #comp_objective = comp_adapt_solution(loc_I_gen, loc_J_gen, agg_supply_bound, demand_bound, cont_perc)
        # comp_objective = 71.2

        c_gen = reshape([norm(loc_I_gen[i,:]-loc_J_gen[j,:]) for j in 1:no_dp for i in 1:no_sp],no_sp,no_dp)

        #static_objective = static_solution(loc_I_gen, loc_J_gen, agg_supply_bound, demand_bound, cont_perc)#


        starting = now()
        o_v, w_v, q_v, p_v, p_true_v = k_adapt_solution(3, loc_I_gen, loc_J_gen, agg_supply_bound, demand_bound, cont_perc)
        ending = now()

        duration = (ending-starting).value/1000 # convert to seconds
        times[n, count] = duration

        #U = enum_uncset(demand_bound, loc_J_gen, cont_perc)
        #Obs = obs_obj(U, find_plan, q_v, c_gen)

        io = open("results/main.txt", "a")
        write(io, ",\t k-adapt overall time: $duration s\n")
        close(io)

    end
end
T = mean(times, dims=2)
plot(1:25, T, ylabel="seconds", xlabel="m", title="Average of 10 computation times for n=5, 3 iterations", label="")
savefig("results/times.png")
#k_curve(o_v, p_v, n_val=no_sp, m_val=no_dp, rel_val=true) 


# k_curve_obs(p_v, Obs, n_val=no_sp, m_val=no_dp)

#obj_f(q, d, c) = 10*sum(max(d[j]-sum(q[i,j] for i in 1:no_sp),0) for j in 1:no_dp) + sum(c[i,j]*q[i,j] for i in 1:no_sp, j in 1:no_dp)
