function random_circuit(g::AbstractGraph, no_cycles::Int)
    k = maximum(degree(g, v) for v in vertices(g))
    ec = edge_color(g, k)
    layers = []
    color_counter = 0

    for i in 1:no_cycles
        color_counter += 1
        color_counter = (color_counter - 1) % k + 1  # 1-based cyclic counter

        layer = []

        # 1-qubit gates: apply a random gate to each site
        append!(layer, (rand(["SqrtX", "SqrtY", "SqrtW"]), [v]) for v in vertices(g))

        # 2-qubit gates: fsim on edges with current color
        for (v1, v2) in ec[color_counter]
            push!(layer, ("fsim", [v1, v2]))
        end

        push!(layers, layer)
    end

    return layers
end
