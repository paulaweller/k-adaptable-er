import Pkg
Pkg.resolve()

include("helpers_data.jl")
include("helpers_solve.jl")

k = parse(Int64, ARGS[1]) 
product = ARGS[2] # supply product
lim = 7200.0         # runtime limit in seconds

read_rio_instance(product)

# # save_x = eval(mod(l, 10) == 0)
save_x = true
results, solutions = run_instance(k, instance, tlim=lim, return_solutions=save_x, pb=false)

# add line/instance number to results
# linedict = Dict(:instance => l)
# merge!(results, linedict)

result_data = DataFrame(results)
#io = open("source/results/rio/individual/results_sol_rio_k$(k).csv", "a")
#close(io)
write_result_to_file("results/rio/$(product)/results_$(product)_k$(k).csv", result_data)

if save_x
     output = open("results/rio/$(product)/individual/sol_rio_$(product)_k$(k).txt", "w")
     CSV.write(output, solutions)
     close(output)
end