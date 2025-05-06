using TensorNetworkQuantumSimulator
const TN = TensorNetworkQuantumSimulator
using ITensorNetworks
const ITN = ITensorNetworks

using NamedGraphs: NamedEdge
using NamedGraphs.GraphsExtensions: add_vertices, add_edges

using Random
using Statistics
using Serialization


function main()
    g = TN.heavy_hexagonal_lattice(2,1)
    g = add_vertices(g, [(0, 3), (3, 4), (7, 4), (10, 1), (7, 0), (3, 0)])
    g = add_edges(g, [NamedEdge((0, 3) => (1, 3)), NamedEdge((3, 4) => (3, 3)), NamedEdge((7, 3) => (7, 4)), NamedEdge((10, 1) => (9, 1)), NamedEdge((7, 0) => (7, 1)), NamedEdge((3, 0) => (3, 1))])

    dict = Dict(
        0 => (0, 3),
        1 => (1, 3),
        2 => (1, 2),
        3 => (1, 1),
        4 => (2, 3),
        5 => (2, 1),
        6 => (3, 4),
        7 => (3, 3),
        8 => (3, 1),
        9 => (3, 0),
        10 => (4, 3),
        11 => (4, 1),
        12 => (5, 3),
        13 => (5, 2),
        14 => (5, 1),
        15 => (6, 3),
        16 => (6, 1),
        17 => (7, 4),
        18 => (7, 3),
        19 => (7, 1),
        20 => (7, 0),
        21 => (8, 3),
        22 => (8, 1),
        23 => (9, 3),
        24 => (9, 2),
        25 => (9, 1),
        26 => (10, 1),
    )
end

main()