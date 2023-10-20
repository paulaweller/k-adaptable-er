using DataFrames, CSV

include("helpers_data.jl")
include("helpers_solve.jl")

sn = 2              # number of service nodes
dn = 4              # number of demand nodes
k = 2               # degree of adaptability
lim = 20.0         # runtime limit in seconds
sed = rand(1:500)   # random seed
inst_gen = generate_instance(sn,dn,sed)
#inst_gen = AllocationInstance([5; 5], [5 10; 2 5], 5, 5, 0.5)

write_instances_to_file([inst_gen], "data/datatest.txt")
#---------------------------------------------------------------
# initialize file

inst_read = read_instances_from_file("data/datatest.txt")

results = run_instance(k, inst_read[1], tlim=lim)

result_data = DataFrame(results)

write_result_to_file("resultdata.csv", result_data)