include("helpers_data.jl")
include("helpers_solve.jl")
include("helpers_plot.jl")

# k=2
# lim = 150.0         # runtime limit in seconds

# inst_read = read_all_instances_from_file("source/data/test_batch.txt")

# results1, solutions1 = run_instance(k, inst_read[1], tlim=lim, return_solutions=true)
# #results2, solutions2 = run_instance(k, inst_read[2], tlim=lim)
# #results3, solutions3 = run_instance(k, inst_read[3], tlim=lim)


# result_data = DataFrame(results1)
# #push!(result_data, results2)
# #push!(result_data, results3)

# write_result_to_file("source/results/resulttest.csv", result_data)
# print(solutions1)
n = 4
m = 15
pc = 0.1
for k = 1:3
    box_plot_from_csv("source/results/data_batch_$(n)_$(m)_$(pc)/combined_results_data_batch_$(n)_$(m)_$(pc)_k$(k)", "data batch $(n)_$(m)_$(pc), k=$k")
end