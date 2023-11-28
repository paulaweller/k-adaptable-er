include("helpers_data.jl")
include("simulate_solution.jl")
include("helpers_plot.jl")

allfile = "source/results/all_batches/combined_results_all_batches_obs"

no = [4,6,8]
mo = [10,15,20]
pco = [0.1]
ko = [2]

# k_plot_from_csv(allfile, method="box", relative=true, observable=true)

# zetas, bbs, boxs = extract_evolutions(no, mo, pco, ko)
# plot_evol(boxs, xlimits=[0,3600], relative=true, name="box_k2_0.1", last = true)

# plot_pc_vs_time(allfile, time=true, objective=true)

# plot_size_vs_time(allfile, percentage=0.1, K=[2])

# plot_zeta_distr_from_csv(allfile, status="unterminated", n_interval=200, relative=false, K=[3,4,5])

# termination_plot_from_csv(allfile, terminated = "bb_infeas", K=[1,2,3,4,5], rel_to_pb=true)


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
