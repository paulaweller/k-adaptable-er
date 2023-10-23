include("helpers_data.jl")
include("helpers_solve.jl")

k=2
lim = 3600.0         # runtime limit in seconds

inst_read = read_instances_from_file("data/test_batch.txt")

results = run_instance(k, inst_read[1], tlim=lim)

result_data = DataFrame(results)

write_result_to_file("results/resultdata.csv", result_data)