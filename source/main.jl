include("helpers_data.jl")
include("simulate_solution.jl")
include("helpers_plot.jl")

allfile = "source/results/all_batches/combined_results_all_batches"

no = [4,6,8]
mo = [10,15,20]
pco = [0.1, 0.3]
ko = [2,3]

# zetas, bbs, boxs = extract_evolutions(no, mo, pco, ko)
# plot_evol(bbs, xlimits=[0,3600], relative=false, name="bb0.3", last = false)

# termination_plot_from_csv(allfile, terminated = "box", k=3, rel_to_pb=true)


for k in ko
    for pc in pco
        println("Solving $(pc)")
        o_pb, o_bb, o_box, oo_pb, oo_box = observable_worst_case_objectives(no, mo, pc, k)
        observables = Dict(
            :o_pb => o_pb, 
            :o_bb => o_bb, 
            :o_box => o_box)
        observable_data = DataFrame(observables)
        write_result_to_file("source/results/all_batches/observables_n$(no)_m$(mo)_pc$(pc)_k$(k).csv", observable_data)

        observable_plot(o_pb,o_bb,o_box,oo_pb, oo_box,k,pc,rel=true)
        observable_plot(o_pb,o_bb,o_box,oo_pb, oo_box,k,pc,rel=false)
    end
end