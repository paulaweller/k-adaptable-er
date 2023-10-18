using LinearAlgebra, Dates, Random, Statistics, JuMP

"""
    AllocationInstance(loc_I,loc_J,W,D,pc)

Instance of the supply preallocation problem where 

    loc_I = locations of service points
    loc_J = locations of demand points
    W     = aggregated supply bound
    D     = maximal demand at one location
    pc    = percentage of affected quantities at risk
"""
struct AllocationInstance
    loc_I::Array{Float64,2}
    loc_J::Array{Float64,2}
    W::Float64
    D::Float64
    pc::Float64
end

"""
    print_instance(instance)

Print instance data.
"""
function print_instance(inst; all=false)
    I = size(inst.loc_I, 2)
    J = size(inst.loc_J, 2)
    s_locations = [(inst.loc_I[1,i], inst.loc_I[2,i]) for i in 1:I]
    d_locations = [(inst.loc_J[1,j], inst.loc_J[2,j]) for j in 1:J]
    println("Instance data:")
    print("\tService locations:")
    for i in 1:I 
        print(" ", s_locations[i])
    end
    print("\n\tDemand locations: ")
    for j in 1:J
        print(" ", d_locations[j])
    end
    print("\n")
    if all==true
        println("\tW = $(inst.W),   D = $(inst.D),   pc = $(inst.pc)")
    end
end


"""
    generate_instance(I_inst, J_inst, seed, demand_bound=5, cont_perc=0.5, agg_supply_bound=round(Int, cont_perc*demand_bound*J), plot_loc=false)

Generate a problem instance with the given parameters. Returns: loc_I_inst, loc_J_inst, demand_bound, cont_perc, agg_supply_bound

Fields:

    I_inst              Number of service points
    J_inst              Number of demand points
    seed                For reproducibility
    demand_bound        Upper bound of demand at one demand point
    cont_perc           Percentage of damage caused by the contingency
    agg_supply_bound    Aggregated supply bound, default: maximal aggregated demand 
    plot_loc            Should the locations of service and demand points be plotted?
    loc_max             How large is the grid
"""
function generate_instance(I_inst, J_inst, seed; demand_bound=5, cont_perc=0.5, agg_supply_bound=round(Int, cont_perc*demand_bound*J_inst), plot_loc=false, loc_max=10)

    # seed for reproducibility
    Random.seed!(seed)
    
    # set of locations for service and demand points
    loc_set = Vector{Float64}[]
    # locations of service points
    loc_I_inst = Vector{Float64}[]
    # locations of demand points
    loc_J_inst = Vector{Float64}[]

    # no two points can have the same location
    while length(loc_set) !== I_inst+J_inst  
        
        # randomly generate locations of service points as x,y-coordinates
        loc_I_inst = rand(1:loc_max,2,I_inst)
        # randomly generate locations of demand points as coordinates
        loc_J_inst = rand(1:loc_max,2,J_inst)
        # add to set (not list) of locations, i.e. a location won't be listed twice
        loc_set = union(union(loc_I_inst[:,i] for i in 1:I_inst),union(loc_J_inst[:,j] for j in 1:J_inst))
    end

    # sort locations by x-coordinate for simplicity
    loc_I_inst = sort(loc_I_inst, dims=2)
    loc_J_inst = sort(loc_J_inst, dims=2)

    instance = AllocationInstance(loc_I_inst,loc_J_inst,agg_supply_bound,demand_bound,cont_perc)
    # plot service and demand points
    if plot_loc == true

        scatter(loc_I_inst[1,:],loc_I_inst[2,:],
                    lims=[0,loc_max+1.2],
                    series_annotations = text.(1:J_inst, :inside),
                    markersize = 20,
                    lab="service point")
        display(scatter!(loc_J_inst[1,:],loc_J_inst[2,:], 
                            series_annotations = text.(1:J_inst, :inside),
                            shape = :square,
                            markersize = 15,
                            lab="demand point",
                            aspect_ratio=:equal
                            ))
    end

    return instance
end

"""
    write_instance_to_file(instances,filename)
"""
function write_instance_to_file(instances::Vector{AllocationInstance}, filename::String)
    # inst = instances[1]
    # no_sp = size(inst.loc_I,2)
    # no_dp = size(inst.loc_J,2)
    # io = open("data/$(filename).txt", "w")
    # write(io, "Pre-Allocation Instances for $(no_sp) service points, $(no_dp) demand points, aggregated demand $(inst.W), maximum demand $(inst.D), affected percentage $(inst.pc)\n")
    # close(io)

    io = open("data/$(filename).txt", "a")
    for instance in instances
        print(io, instance, "\n")    
    end
    close(io)
end

"""
    read_instance_from_file(filename)

returns a vector with instances
"""
function read_instance_from_file(filename::String)

    instances = AllocationInstance[]
    io = open(filename)

    for ln in eachline(io)
        inst = eval(Meta.parse(ln))
        push!(instances, inst)
    end
    close(io)

    return instances
end