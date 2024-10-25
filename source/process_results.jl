using DataFrames, CSV

include("simulate_solution.jl")

# # process synthetic data results 

# extract_results([4,6,8],[10,15,20], [0.1, 0.3], [1,2,3,4,5])
# observ = add_observables([4,6,8], [10,15,20], [0.1,0.3], [1])
# results = DataFrame(CSV.File("source/results/all_batches/combined_results_all_batches.csv"))
# results = innerjoin(results, observ, on = [:n,:m,:k,:pc,:instance])
# output = open("source/results/all_batches/combined_results_all_batches_obs_k1.csv", "w")
# CSV.write(output, results)
# close(output)

# process rio results

ko = [1,2,3,4,5]
supplies = ["Food", "Water", "Hygiene", "Cleaning", "Mattress", "Medicine"]
extract_results(ko, supplies)
observ = add_observables(supplies, ko)
results = DataFrame(CSV.File("source/results/rio/combined_results.csv"))
results = innerjoin(results, observ, on = [:supply,:k])
output = open("source/results/rio/combined_results_rio_obs.csv", "w")
CSV.write(output, results)
close(output)

