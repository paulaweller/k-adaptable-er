using Gurobi, JuMP

model = Model(Gurobi.Optimizer)

@variable(model, x[1:1,1:3], Bin)

i=length(x)+1
push!(x, @variable(model, base_name = "x[$i,$k]", binary=true))

@objective(model, Max, x[1]+2x[2])




print(model)