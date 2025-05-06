using NamedGraphs: NamedGraph
include("utils.jl")

function heavyhex_isingspinglass_without_zzz_circuit_layer(g::NamedGraph, γ::Float64, β::Float64, dv, dij)
    ec = edge_color(g, 3)
    layers = []
    layer = []
    # For 1-qubit gates
    for v in vertices(g)
        push!(layer, ("Rz", [v], 2*γ*dv[v]))
    end

    # For 2-qubit gates
    for c=1:3
        edgegroup = ec[c]
        for e in edgegroup
            srce, dste = first(e), last(e)
            es = NamedEdge(srce => dste)
            if es in keys(dij)
                push!(layer, ("Rzz", [src(es), dst(es)], 2*γ*dij[es]))
            end
        end
    end

    for v in vertices(g)
        push!(layer, ("Rx", [v], 2*β))
    end

    push!(layers, layer)
    return layers
end
