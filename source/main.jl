include("helpers_data.jl")
#include("helpers_solve.jl")
include("helpers_plot.jl")
using JSON

allfile = "source/results/all_batches/combined_results_all_batches"

no = [4,6,8]
mo = [10,15,20]
pco = [0.1, 0.3]
ko = [2]
# zetas, bbs, boxs = extract_evolutions(no, mo, pco, ko)

for pc in pco
    plot_size_vs_time(allfile, percentage=pc, K=[3])
end