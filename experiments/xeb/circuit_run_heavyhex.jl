using TensorNetworkQuantumSimulator
const TN = TensorNetworkQuantumSimulator
using ITensorNetworks
const ITN = ITensorNetworks

using NamedGraphs: NamedEdge
using NamedGraphs.GraphsExtensions: add_vertices, add_edges

using Random
using Statistics
using Serialization
Random.seed!(1734)

include("circuit_generators.jl")

function main()
    g = TN.heavy_hexagonal_lattice(2,1)
    g = add_vertices(g, [(0, 3), (3, 4), (7, 4), (10, 1), (7, 0), (3, 0)])
    g = add_edges(g, [NamedEdge((0, 3) => (1, 3)), NamedEdge((3, 4) => (3, 3)), NamedEdge((7, 3) => (7, 4)), NamedEdge((10, 1) => (9, 1)), NamedEdge((7, 0) => (7, 1)), NamedEdge((3, 0) => (3, 1))])

    sites = ITN.siteinds("S=1/2", g)

    ψ0 = ITN.ITensorNetwork(v -> "↑", sites)

    maxdim, cutoff = 64, 1e-14
    apply_kwargs = (; maxdim, cutoff, normalize = true)
    #Parameters for BP, as the graph is not a tree (it has loops), we need to specify these
    set_global_bp_update_kwargs!(;
        maxiter = 30,
        tol = 1e-10,
        message_update_kwargs = (;
            message_update_function = ms -> make_eigs_real.(ITN.default_message_update(ms))
        ),
    )

    ψ = copy(ψ0)
    ψψ = build_bp_cache(ψ)

    circuit = random_circuit(g, 5)

    for (cycle, layer) in enumerate(circuit)
        ψ, ψψ, errors = apply(layer, ψ, ψψ; apply_kwargs, verbose = false)
    end
end

main()