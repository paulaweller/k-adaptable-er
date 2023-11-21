using DataFrames, CSV

function extract_results(no, mo, pco, ko)
    results = DataFrame()
    for n in no
        for m in mo
            for pc in pco
                for k in ko
                    # Specify the directory path
                    directory_path = "source/results/data_batch_$(n)_$(m)_$(pc)"

                    # Specify the common beginning of the filename
                    common_prefix = "combined_results_data_batch_$(n)_$(m)_$(pc)_k$(k).csv"

                    # Get a list of files in the directory that match the pattern
                    files = filter(x -> occursin(common_prefix, x), readdir(directory_path))
                    for file in files
                        # get file content
                        data = DataFrame(CSV.File(directory_path*"/"*file))

                        # add columns for n, m, pc, k
                        data[!,:n] .= n
                        data[!,:m] .= m
                        data[!,:pc] .= pc
                        data[!, :k] .= k

                        #merge with other data
                        results = vcat(results, data)
                    end
                end
            end
        end
    end
    output = open("source/results/all_batches/combined_results_all_batches.csv", "w")
    CSV.write(output, results)
    close(output)
    return results
end

# extract_results([4,6,8],[10,15,20], [0.1, 0.3], [1,2,3,4,5])

function add_observables(no, mo, pco, ko)
    
    observ = DataFrame()

    for pc in pco
        println(" pc = $(pc)")
        for k in ko
            println("\n k = $k")
            obs_df = observable_worst_case_objectives(no, mo, pc, k)

            #merge with other data
            observ= vcat(observ, obs_df)
        end
    end

    
    return observ
end

observ = add_observables([4,6,8], [10,15,20], [0.1,0.3], [1,2,3,4,5])
results = DataFrame(CSV.File("source/results/all_batches/combined_results_all_batches.csv"))
results = innerjoin(results, observ, on = [:n,:m,:k,:pc,:instance])
    output = open("source/results/all_batches/combined_results_all_batches_obs.csv", "w")
    CSV.write(output, results)
    close(output)