include("helpers_data.jl")
include("helpers_solve.jl")

k=2
lim = 150.0         # runtime limit in seconds

inst_read = read_instances_from_file("data/test_batch.txt")

results1, solutions1 = run_instance(k, inst_read[1], tlim=lim)
#results2, solutions2 = run_instance(k, inst_read[2], tlim=lim)
#results3, solutions3 = run_instance(k, inst_read[3], tlim=lim)


result_data = DataFrame(results1)
#push!(result_data, results2)
#push!(result_data, results3)

write_result_to_file("results/resultdatabatchnewnew.csv", result_data)