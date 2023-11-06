include("helpers_data.jl")
include("helpers_solve.jl")

k = parse(Int64, ARGS[1]) 
l = parse(Int64, ARGS[2]) #task id
file = ARGS[3] #data file
lim = 150.0         # runtime limit in seconds

inst_read = read_one_instance_from_file(l, "data/$(file).txt")

# TODO TEST THIS
save_x = eval(mod(l, 5) == 0)
results, solutions = run_instance(k, inst_read, tlim=lim, return_solutions=save_x)

#results, solutions = run_instance(k, inst_read, tlim=lim)

result_data = DataFrame(results)
write_result_to_file("results/$(file)/results_$(file)_k$(k)_l$(l).csv", result_data)

# TODO TEST THIS
if save_x
     output = open("results/$(file)/sol_$(file)_k$(k)_l$(l).txt", "w")
     CSV.write(output, solutions)
     close(output)
 end