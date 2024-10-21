import Pkg
Pkg.resolve()

include("helpers_data.jl")
include("helpers_solve.jl")

k = parse(Int64, ARGS[1]) 
product = parse(String, ARGS[2]) # supply product
lim = 7200.0         # runtime limit in seconds

# supplies
W_df = CSV.read("data/rio_supplies.csv", DataFrame)
W = W_df[W_df[!,1].=="persons total", product][1]

# demand bound per location
D_df = CSV.read("data/rio_demand.csv", DataFrame)
max_pc_per_area = 0.3679
D = D_df[D_df[!,1].!="Total", "Population"].*max_pc_per_area

# percentage for aggregate demand bound
pc_overall = 0.1096

# distances between supply and demand points
dist_df = CSV.read("data/rio_dist.csv", DataFrame)
dist = Array(dist_df[!,2:end])

# empty set of locations
loc = Array{Float64}(undef, 0, 0)

instance = AllocationInstance(loc,loc,W,D,pc_overall,dist)

# # save_x = eval(mod(l, 10) == 0)
save_x = true
results, solutions = run_instance(k, instance, tlim=lim, return_solutions=save_x)

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