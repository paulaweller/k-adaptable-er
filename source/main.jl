include("helpers_data.jl")
#include("helpers_solve.jl")
include("helpers_plot.jl")
using JSON

allfile = "source/results/all_batches/combined_results_all_batches"

no = [4,6,8]
mo = [10,15,20]
pco = [0.1, 0.3]
ko = [1,2,3]
#zetas, bbs, boxs = extract_evolutions(no, mo, pco, ko)


#plot_evol(boxs, xlimits=[0,600], relative=true, name="box", last = true)

termination_plot_from_csv(allfile, terminated = "box", k=3, rel_to_pb=true)