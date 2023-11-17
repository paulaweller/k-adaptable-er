include("helpers_data.jl")
include("simulate_solution.jl")
include("helpers_plot.jl")
using JSON

allfile = "source/results/all_batches/combined_results_all_batches"

no = [4,6,8]
mo = [10,15,20]
pco = [0.1,0.3]
ko = [1,2,3]
#zetas, bbs, boxs = extract_evolutions(no, mo, pco, ko)


#plot_evol(boxs, xlimits=[0,600], relative=true, name="box", last = true)

# termination_plot_from_csv(allfile, terminated = "box", k=3, rel_to_pb=true)



function observable_plot(o_pb, o_bb, o_box, oo_pb, oo_box,k, pc; rel=true)
    if rel == true
        plot1 = plot(ylabel="observable improvement")
        boxplot!(plot1, ["bb ($(length(o_bb)))" "box"], [100 .*(o_pb.-o_bb)./o_pb 100 .*(o_pb.-o_box)./o_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        if length(oo_pb) > 0
            boxplot!(plot1, ["bb infeas ($(length(oo_box)))"], [100 .*(oo_pb.-oo_box)./oo_pb], 
                    leg=false, 
                    linewidth=2,
                    linecolour= :match,
                    fillalpha = 0.4)
        end
        savefig(plot1, "source/plots/observable/rel_objective_$(pc)_k$(k).pdf")
    elseif rel == false
        plot1 = plot(ylabel="observable worst-case objective")
        boxplot!(plot1, ["bb($(length(o_bb)))" "box" "pb"], [o_bb o_box o_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        if length(oo_pb) > 0
                boxplot!(plot1, ["bb infeas ($(length(oo_box))): box" " bb infeas: pb"], [oo_box oo_pb], 
                leg=false, 
                linewidth=2,
                linecolour= :match,
                fillalpha = 0.4)
        end
        savefig(plot1, "source/plots/observable/objective_$(pc)_k$(k).pdf")
    else error("optinal argument rel must be true or false")
    end
end


for k in ko
    for pc in pco
        o_pb, o_bb, o_box, oo_pb, oo_box = observable_worst_case_objectives(no, mo, pc, k)
        observable_plot(o_pb,o_bb,o_box,oo_pb, oo_box,k,pc,rel=false)
    end
end