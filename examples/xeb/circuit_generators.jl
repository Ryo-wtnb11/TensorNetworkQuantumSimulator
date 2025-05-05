function random_circuit(g::AbstractGraph, no_cycles::Int64)
    k = maximum([degree(g, v) for v in vertices(g)])
    ec = edge_color(g, k)
    gates = []
    color_counter = 0
    for i = 1:no_cycles
        color_counter += 1
        single_site_gates = [(rand(["SqrtX", "SqrtY", "SqrtW"]), v) for v in vertices(g)]
        two_site_gates = [
            ("fsim", src(first(e)), dst(first(e))) for
            e in filter(e -> last(e) == color_counter, ec)
        ]
        push!(gates, (single_site_gates, two_site_gates))
        color_counter = color_counter % k
    end
    return gates
end

function SimpleGraphAlgorithms.edge_color(g::NamedGraph, k::Int64)
    pg, vs = position_graph(g), collect(vertices(g))
    ec = edge_color(UG(pg), k)
    ec_g = [NamedEdge(vs[first(first(e))] => vs[last(first(e))]) => last(e) for e in ec]
    return ec_g
end
