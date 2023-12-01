include("helpers_data.jl")
include("simulate_solution.jl")
include("helpers_plot.jl")

allfile = "source/results/all_batches/combined_results_all_batches_obs"

k_plot_from_csv(allfile, method="box", relative=false, observable=false)

# zetas, bbs, boxs = extract_evolutions([4,6,8], [10,15,20], [0.1,0.3], [2])
# plot_evol(boxs, allfile, xlimits=[0,3600], relative=true, name="box_k2", last = true)

# plot_pc_vs_time(allfile, time=true, objective=true)

#plot_size_vs_time(allfile, percentage=0.1, K=[1,2])

# plot_zeta_distr_from_csv(allfile, status="terminated", n_interval=300, relative=false, K=[1,2,3])

# termination_plot_from_csv(allfile, terminated = "neither", K=[1,2,3,4,5], rel_to_pb=false)


# for k in ko
#     for pc in pco
#         println("Solving $(pc)")
#         o_pb, o_bb, o_box, oo_pb, oo_box = observable_worst_case_objectives(no, mo, pc, k)
#         observables = Dict(
#             :o_pb => o_pb, 
#             :o_bb => o_bb, 
#             :o_box => o_box)
#         observable_data = DataFrame(observables)
#         write_result_to_file("source/results/all_batches/observables_n$(no)_m$(mo)_pc$(pc)_k$(k).csv", observable_data)

#         observable_plot(o_pb,o_bb,o_box,oo_pb, oo_box,k,pc,rel=true)
#         observable_plot(o_pb,o_bb,o_box,oo_pb, oo_box,k,pc,rel=false)
#     end
# end
