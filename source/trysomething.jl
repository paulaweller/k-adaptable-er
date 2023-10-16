include("helpers_data.jl")

i = [generate_instance(2,3,i) for i in 1:2]

write_instance_to_file(i, "datatest")

ins = read_instance_from_file("data/datatest.txt")

for ini in ins
    print_instance(ini, all=true)
end
