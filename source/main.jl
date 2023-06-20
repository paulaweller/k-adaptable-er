using LinearAlgebra, StatsPlots

include("helpers.jl")
#include("solve_iter.jl")
include("solve_bb.jl")
include("solve_boxidea_inplace.jl")
include("solve_bb_inplace.jl")
include("solve_boxidea.jl")


# files = ["results/larger_inst25.txt", "results/larger_inst50.txt", "results/larger_inst75.txt", "results/larger_inst100.txt"]
# lbls = ["25" "50" "75" "100"]
# box_plot_from_files(files, lbls)

# seed 9,7 is great
#61 weird
@show inst_gen = generate_instance(1,3,22)
#inst_gen = AllocationInstance([2 4; 8 3; 9 7], [1 9; 3 2; 4 1; 6 3; 7 9; 9 10], 15, 5, 0.5)
#inst_gen = AllocationInstance([5 5], [5 2; 10 5], 5, 5, 0.5)

I = size(inst_gen.loc_I, 1)
J = size(inst_gen.loc_J, 1)
c = reshape([norm(inst_gen.loc_I[i,:]-inst_gen.loc_J[j,:]) for j in 1:J for i in 1:I],I,J)
# println("Locations:\nloc_I = $(inst_gen.loc_I)\nloc_J = $(inst_gen.loc_J)")
#----------------------------

#theta_iter, x_iter, y_iter, p_iter, p_true_iter = k_adapt_solution(2, inst_gen)
# starting = now()

######## x_general, y_general, s_general, theta_general, it_general, runtime_general = solve_bb_general(2, inst_gen)

# ending = now()
# duration_general = (ending-starting).value/1000 # convert to seconds

######## x_gen2, y_gen2, s_gen2, theta_gen2, it_gen2, runtime_gen2 = solve_bb_inplace(2, inst_gen)

@show x_box, y_box, s_box, xi_box, theta_box, it_box, runtime_box = solve_boxes(2, inst_gen)

x_box2, y_box2, s_box2, xi_box2, theta_box2, it_box2, runtime_box2 = solve_boxes_inplace(2, inst_gen)

# println("instance: $(inst_gen)")

# println("x_general = $(x_general)")
# println("y_general = $(y_general)")
# println("s_general = $(s_general)")

# println("x_box = $(x_box)")
# println("y_box = $(y_box)")
# println("s_box = $(s_box)")

# #obj_val, w_val, q_val, p_val, p_true =  k_adapt_solution(7, inst_gen)

# x_diff = x_general - x_box
# y_diff = y_general - y_box
# println("BB vs. box idea solution:\n 
#             Objective: $(theta_general)   |    $(theta_box)\n
#             Time: $(duration_general)   |    $(duration_box)\n
#             Iterations: $(it_general)   |    $(it_box)\n
#             Difference in x: $(x_diff)\n
#             Difference in y: $(y_diff) ")
#----------------------------
# obj_val = theta_iter[end]
# w_val = x_iter[end]
# q_val = y_iter[end]

# println("k-adaptable solution")
# println("objectives: $obj_val")
# println("supply status = $w_val")
# println("transportation: $(q_val)")

#k_curve(obj_val, p_val, comp_adapt = 71.2)
#observations = obs_obj(U, find_plan, q_val, c)
# k_curve_obs(p_val, observations) #with W = 1