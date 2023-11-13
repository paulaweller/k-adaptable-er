include("helpers_data.jl")
include("helpers_solve.jl")
include("helpers_plot.jl")

# k=2
# lim = 150.0         # runtime limit in seconds

# inst_read = read_all_instances_from_file("source/data/test_batch.txt")

# results1, solutions1 = run_instance(k, inst_read[1], tlim=lim, return_solutions=true)
# #results2, solutions2 = run_instance(k, inst_read[2], tlim=lim)
# #results3, solutions3 = run_instance(k, inst_read[3], tlim=lim)


# result_data = DataFrame(results1)
# #push!(result_data, results2)
# #push!(result_data, results3)

# write_result_to_file("source/results/resulttest.csv", result_data)
# print(solutions1)


no = [8]
mo = [10,15,20]
pco = [0.1, 0.3]
ko = [1,2,3]
for k in ko
for n in no
    for m in mo
        for pc in pco
            # box_plot_from_csv("source/results/data_batch_$(n)_$(m)_$(pc)/combined_results_data_batch_$(n)_$(m)_$(pc)_k$(k)", "data batch $(n)_$(m)_$(pc), k=$k", times=false)
            filename = "source/results/all_batches/combined_results_data_batch_$(n)_$(m)_$(pc)_k$(k)"
            data = DataFrame(CSV.File("$(filename).csv"))
            data[!, :n] .= n
            data[!, :m] .= m
            data[!, :pc] .= pc
            data[!, :k] .= k
            output = open("$(filename).csv", "w")
            CSV.write(output, data)
            close(output)
        end
    end
end
end

