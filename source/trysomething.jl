include("helpers_data.jl")

for l in 1:3
    inst_read = read_one_instance_from_file(l, "source/data/test_batch.txt")
    print_instance(inst_read)
end

