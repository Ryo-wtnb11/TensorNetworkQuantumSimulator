using TensorNetworkQuantumSimulator
const TN = TensorNetworkQuantumSimulator
using ITensorNetworks
const ITN = ITensorNetworks

using Random
using Statistics
using Serialization
Random.seed!(1734)

include("circuit_generators.jl")

function main()
    g = sycamore_53_qubit_grid(;patch=1)

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
    @show circuit
    return 0

    for (cycle, layer) in enumerate(circuit)
        ψ, ψψ, errors = apply(layer, ψ, ψψ; apply_kwargs, verbose = false)
        serialize("ψ_$(cycle).jld", ψ)
        serialize("errors_$(cycle).jld", errors)
    end
end

main()