include("helpers_data.jl")
include("helpers_solve.jl")

lim = 240.0         # runtime limit in seconds

inst_read = read_instances_from_file("data/datatest.txt")

results = run_instance(k, inst_read[1], tlim=lim)

result_data = DataFrame(results)

write_result_to_file("results/resultdata_triton.csv", result_data)