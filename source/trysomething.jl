
 N = []
 tau = [[],[],[]]
 xi = [1,1]
for k in 1:3
    # each child node is the current uncset configuration...
    tau_temp = copy(tau)
    # ...with the new scenario added to uncset number k
    # tau_temp[k] = union(tau_temp[k], [xi])
    append!(tau_temp[k], [xi])
    println(tau_temp)
    #N = union(N, [tau_temp])
    @show push!(N, tau_temp) 
end