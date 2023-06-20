using Gurobi, JuMP

function create_model()
    model = Model(Gurobi.Optimizer)

    #@variable(model, x[1:2,1:3])

    x = model[:x] = @variable(model, [1:1,1:3], base_name = "x")

    @constraint(model, demi[i=1:2], x[i]+x[i+1] <= 10)

    @objective(model, Max, x[1]+2x[2])

    return model
end

function modify_model(model)

    t = size(model[:x],1)+1
    x_new = @variable(model, [[t],1:3], base_name = "x")
    model[:x] = [model[:x]; x_new]

    return model[:x]
end

modeli = create_model()

modify_model(modeli)
modify_model(modeli)