include("helpers_data.jl")
#include("helpers_solve.jl")
include("helpers_plot.jl")
using JSON

allfile = "source/results/all_batches/combined_results_all_batches"


function extract_zetas()
    no = [4,6,8]
    mo = [10,15,20]
    pco = [0.1, 0.3]
    ko = [1,2,3]
    zetas = Dict()
    for n in no
        for m in mo
            for pc in pco
                for k in ko
                    # Specify the directory path
                    directory_path = "source/results/data_batch_$(n)_$(m)_$(pc)/individual"

                    # Specify the common beginning of the filename
                    common_prefix = "sol_data_batch_$(n)_$(m)_$(pc)_k$(k)_"

                    # Get a list of files in the directory that match the pattern
                    files = filter(x -> occursin(common_prefix, x), readdir(directory_path))
                    i = 0
                    for file in files
                        i = i+1
                        println(file)
                        dict = read_solutions_from_file("source/results/data_batch_$(n)_$(m)_$(pc)/individual/$(file)")
                        dict_zeta = Dict([n,m,pc,k,i] => dict["evol_zeta"])
                        merge!(zetas, dict_zeta)
                    end
                end
            end
        end
    end
    return zetas
end

# Initialize an empty plot
function vecfromstr(str)
    str = split(chop(str, head=1), ",")
    str = strip.(str, [' '])
    str = strip.(str, ['['])
    str = strip.(str, [']'])
    str_s = parse.(Float64,str)
    final = transpose(reshape(str_s, 2, :))
    return final
end

# Iterate over dictionary keys and values
function plot_z(zetas)

    time_plot = plot()
    for (key, values) in zetas
        # Extract x and y coordinates from each 2D vector

        coordinates = vecfromstr(values)

        # Extract x and y coordinates from each 2D vector
        x_vals = coordinates[:,1]
        y_vals = coordinates[:,2]
        # Plot the data
        plot!(time_plot, x_vals, y_vals)
    end
    savefig(time_plot, "source/plots/zetas/zeta.pdf")

end

zz = extract_zetas()
plot_z(zz)