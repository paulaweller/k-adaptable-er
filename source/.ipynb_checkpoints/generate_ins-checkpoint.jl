include("helpers_data.jl")

sn_set = [4,6,8]           # number of service nodes
dn_set = [10, 15,20]       # number of demand nodes
pc_set = [0.1, 0.3]
i_n = 50                   # number of instances
# initialize seed file
io = open("data/seeds.txt", "w")
close(io)

for sn in sn_set
    for dn in dn_set
        for pc in pc_set
            #random seed
            sed      = rand(1:10000, i_n)
            # save seeds
            io = open("data/seeds.txt", "a")
            print(io, "sn=$(sn), dn=$(dn), pc=$(pc)", sed, "\n")
            close(io)
            # initialize file with empty instance
            filename = "data/data_batch_$(sn)_$(dn)_$(pc).txt"
            write_instances_to_file(AllocationInstance[], filename)
            # generate and write instances to file
            for i in 1:i_n
                inst_gen = generate_instance(sn,dn,sed[i],cont_perc=pc)
                add_instances_to_file([inst_gen], filename)
            end
        end
    end
end
