using LinearAlgebra #, DataFrames

include("helpers_data.jl")
include("helpers_solve.jl")

sn = 2              # number of service nodes
dn = 5              # number of demand nodes
k = 3               # degree of adaptability
lim = 250.0         # runtime limit in seconds
sed = rand(1:500)   # random seed
inst_gen = generate_instance(sn,dn, sed)
#inst_gen = AllocationInstance([5; 5], [5 10; 2 5], 5, 5, 0.5)

#---------------------------------------------------------------

run_instance(k, inst_gen, tlim=120.0)

nothing