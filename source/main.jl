include("helpers_data.jl")
include("simulate_solution.jl")
include("helpers_plot.jl")
using JSON

allfile = "source/results/all_batches/combined_results_all_batches"

no = [4,6,8]
mo = [10,15,20]
pco = [0.1]
ko = [4]

# zetas, bbs, boxs = extract_evolutions(no, mo, pco, ko)
# plot_evol(bbs, xlimits=[0,3600], relative=false, name="bb0.3", last = false)

# termination_plot_from_csv(allfile, terminated = "box", k=3, rel_to_pb=true)


for k in ko
    for pc in pco
        o_pb, o_bb, o_box, oo_pb, oo_box = observable_worst_case_objectives(no, mo, pc, k)
        observable_plot(o_pb,o_bb,o_box,oo_pb, oo_box,k,pc,rel=false)
    end
end