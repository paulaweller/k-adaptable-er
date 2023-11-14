include("helpers_data.jl")
#include("helpers_solve.jl")
include("helpers_plot.jl")

# plotting times and objectives
# opt = ["both", "bb", "box", "neither", "bb_infeas"]
# K = [2,3]
# for ko in K
#     for o in opt
#         k_plot_from_csv("source/results/all_batches/combined_results_all_batches", method="box")
#     end
# end

k_plot_from_csv("source/results/all_batches/combined_results_all_batches", method="bb", relative=false)