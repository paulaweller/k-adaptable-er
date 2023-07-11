using JuMP

function create_model()
    model = Model()
end

function modify_model!(model)
    z = @variable(model, [1:2])
    @constraint(model, z[1]+z[2] <= 4)
end

m = create_model()
print(m)
modify_model!(m)
modify_model!(m)
print(m)