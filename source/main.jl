using LinearAlgebra, DataFrames, CSV

include("helpers_data.jl")
include("helpers_solve.jl")

sn = 1              # number of service nodes
dn = 2              # number of demand nodes
k = 2               # degree of adaptability
lim = 250.0         # runtime limit in seconds
sed = rand(1:500)   # random seed
inst_gen = generate_instance(sn,dn, sed)
#inst_gen = AllocationInstance([5; 5], [5 10; 2 5], 5, 5, 0.5)

#---------------------------------------------------------------
# initialize file

results = run_instance(k, inst_gen, tlim=180.0, bb=false)

result_data = DataFrame(results)

output = open("results/resultdata.csv", "a")
CSV.write(output, result_data)
close(output)

output = open("results/resultdata.csv", "a")
CSV.write(output, result_data, append=true)
close(output)