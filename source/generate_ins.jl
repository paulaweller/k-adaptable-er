include("helpers_data.jl")

sn = 4              # number of service nodes
dn = 8              # number of demand nodes

sed = rand(1:500)   # random seed
inst_gen = generate_instance(sn,dn,sed)
#inst_gen = AllocationInstance([5; 5], [5 10; 2 5], 5, 5, 0.5)

write_instances_to_file([inst_gen], "data/datatest.txt")