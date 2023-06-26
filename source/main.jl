using LinearAlgebra, ProfileView

include("helpers.jl")
include("solve_iter.jl")
include("solve_bb.jl")
include("solve_boxes_inplace.jl")
include("solve_bb_inplace.jl")
include("solve_boxes.jl")


# files = ["results/larger_inst25.txt", "results/larger_inst50.txt", "results/larger_inst75.txt", "results/larger_inst100.txt"]
# lbls = ["25" "50" "75" "100"]
# box_plot_from_files(files, lbls)

# seed 9,7 is great
#61 weird
sn = 1
dn = 3
k = 3
@show inst_gen = generate_instance(1,3,850)
#inst_gen = AllocationInstance([5; 5], [5 10; 2 5], 5, 5, 0.5)
it = round(Int, 1 + (k-1)/(dn-1))

#----------------------------
# @code_warntype k_adapt_solution(it, inst_gen)

# @show theta_iter, x_iter, y_iter, p_iter, p_true_iter = k_adapt_solution(it, inst_gen)

@show x_box, y_box, s_box, xi_box, theta_box, it_box, runtime_box = solve_boxes(k, inst_gen)

@show x_box2, y_box2, s_box2, xi_box2, theta_box2, it_box2, runtime_box2 = solve_boxes_inplace(k, inst_gen)

@show x_general, y_general, s_general, theta_general, it_general, runtime_general = solve_bb_general(k, inst_gen)

@show x_gen2, y_gen2, s_gen2, theta_gen2, it_gen2, runtime_gen2 = solve_bb_inplace(k, inst_gen)

nothing
# x_diff = x_general - x_box
# y_diff = y_general - y_box
# println("BB vs. box idea solution:\n 
#             Objective: $(theta_general)   |    $(theta_box)\n
#             Time: $(duration_general)   |    $(duration_box)\n
#             Iterations: $(it_general)   |    $(it_box)\n
#             Difference in x: $(x_diff)\n
#             Difference in y: $(y_diff) ")
#----------------------------