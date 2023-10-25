include("helpers_data.jl")

sn = 4              # number of service nodes
dn = 8              # number of demand nodes
i_n = 3              # number of instances

sed = rand(1:500, i_n)   # random seed
filename = "data/test_batch.txt"

# initialize file 
write_instances_to_file(AllocationInstance[], filename)

for i in 1:i_n
    inst_gen = generate_instance(sn,dn,sed[i])
    add_instances_to_file([inst_gen], filename)
end
