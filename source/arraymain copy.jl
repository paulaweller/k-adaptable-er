include("helpers_data.jl")
include("helpers_solve.jl")

k = 2
l = 6 #task id
file = "source/data/data_batch_4_10_0.3.txt" #data file
lim = 3600.0         # runtime limit in seconds

inst_read = read_one_instance_from_file(l, file)

# TODO TEST THIS
# save_x = eval(mod(l, 10) == 0)
save_x = false
results, solutions = run_instance(k, inst_read, tlim=lim, return_solutions=save_x,bb=false,pb=false)

# add line/instance number to results
linedict = Dict(:instance => l)
merge!(results, linedict)

result_data = DataFrame(results)
write_result_to_file("source/results/nosym/results_k$(k)_l$(l)_sym.csv", result_data)

# TODO TEST THIS
if save_x
     output = open("results/nosym/$(file)/individual/sol_$(file)_k$(k)_l$(l).txt", "w")
     CSV.write(output, solutions)
     close(output)
end