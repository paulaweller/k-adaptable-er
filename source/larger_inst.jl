using LinearAlgebra, Plots, PrettyTables, Statistics

include("solve_iter.jl")
include("helpers.jl")
no_sp = 5
no_dp = 30
no_sample = 10
times = zeros(no_sample)
values = zeros(no_sample)
io = open("results/larger_inst$no_dp.txt", "w")
write(io, "K-adapt times vs objective % for $no_sp service points, $no_dp demand points, 2 iterations\n")
close(io)
rand_seeds = rand(500:999, no_sample)

for s in 1:no_sample

    se = rand_seeds[s]
    io = open("results/larger_inst$no_dp.txt", "a")
    write(io, "service points: $no_sp,\t demand points: $no_dp,\t seed: $se")
    close(io)
    #se=8

    inst = generate_instance(no_sp,no_dp, se, cont_perc=0.1, plot_loc=false, loc_max=20)
    #write(io, "loc_I_gen = $loc_I_gen, loc_J_gen = $loc_J_gen, demand_bound = $demand_bound, cont_perc = $cont_perc \n")


    starting = now()
    o_v, w_v, q_v, p_v, p_true_v = k_adapt_solution(2, inst)
    ending = now()

    duration = (ending-starting).value/1000 # convert to seconds
    times[s] = duration
    values[s] = 100*o_v[2]/o_v[1]

    #U = enum_uncset(demand_bound, loc_J_gen, cont_perc)
    #Obs = obs_obj(U, find_plan, q_v, c_gen)

    io = open("results/larger_inst$no_dp.txt", "a")
    write(io, ",\t k-adapt overall time: $duration s, objective % $(values[s])\n")
    close(io)

end


scatter(times, values, ylabel="objective %", xlabel="time in seconds", title="K-adaptable solutions after 2 iterations, n=$no_sp, m=$no_dp, p=0.1", label="")
savefig("results/larger_inst$no_dp.pdf")
#k_curve(o_v, p_v, n_val=no_sp, m_val=no_dp, rel_val=true) 


# k_curve_obs(p_v, Obs, n_val=no_sp, m_val=no_dp)

#obj_f(q, d, c) = 10*sum(max(d[j]-sum(q[i,j] for i in 1:no_sp),0) for j in 1:no_dp) + sum(c[i,j]*q[i,j] for i in 1:no_sp, j in 1:no_dp)
