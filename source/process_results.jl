using DataFrames, CSV

include("simulate_solution.jl")

function extract_results_rio(ko,supplies)
    results = DataFrame()

    for supply in supplies

        # Specify the directory path
        directory_path = "source/results/rio/$(supply)"

        for k in ko
            # Specify the common beginning of the filename
            common_prefix = "results_$(supply)_k$(k).csv"

            # Get a list of files in the directory that match the pattern
            files = filter(x -> occursin(common_prefix, x), readdir(directory_path))
            for file in files
                # get file content
                data = DataFrame(CSV.File(directory_path*"/"*file))

                # add columns for supply, k
                data[!,:supply] .= supply
                data[!, :k] .= k

                #merge with other data
                results = vcat(results, data)
            end
        end
    end

    output = open("source/results/rio/combined_results.csv", "w")
    CSV.write(output, results)
    close(output)
    return results
end

function add_observables_rio(supplies, ko)
    
    observ = DataFrame()

    for supply in supplies
        println(" supply = $(supply)")
        for k in ko
            println("\n k = $k")
            obs_df = observable_worst_case_objectives(k,supply)

            #merge with other data
            observ= vcat(observ, obs_df)
        end
    end

    
    return observ
end

ko = [1,2,3,4,5]
supplies = ["Food", "Water", "Hygiene", "Cleaning", "Mattress", "Medicine"]
extract_results(ko, supplies)
observ = add_observables(supplies, ko)
results = DataFrame(CSV.File("source/results/rio/combined_results.csv"))
results = innerjoin(results, observ, on = [:supply,:k])
output = open("source/results/rio/combined_results_rio_obs.csv", "w")
CSV.write(output, results)
close(output)

