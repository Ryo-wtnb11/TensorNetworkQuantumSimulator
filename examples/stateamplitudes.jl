using Random
using Statistics
using Dictionaries

using TensorNetworkQuantumSimulator
const TN = TensorNetworkQuantumSimulator

using ITensorNetworks: ITensorNetworks, siteinds, Algorithm
const ITN = ITensorNetworks
using NamedGraphs.NamedGraphGenerators: named_grid
using Statistics

using ITensors

using StatsBase

Random.seed!(1734)

using Serialization

function main()
    g = TN.heavy_hexagonal_lattice(1,1)
    s = siteinds("S=1/2", g)
    ψ = ITN.random_tensornetwork(ComplexF64, s; link_space = 2)
    ψ, ψψ = normalize(ψ)

    n = 1 # You can change this to sample multiple amplitudes

    set_global_bp_update_kwargs!(;
        maxiter = 30,
        tol = 1e-10,
        message_update_kwargs = (;
            message_update_function = ms -> make_eigs_real.(ITN.default_message_update(ms))
        ),
    )

    set_global_boundarymps_update_kwargs!(
        message_update_kwargs = (; niters = 25, tolerance = 1e-12),
    )

    probabilities_bp = []
    probabilities_bmps = []
    probabilities_loop = []

    for sample=1:n
        bit_string = Dictionary(vertices(ψ), [rand([1, 2]) for v in vertices(ψ)])
        pψ = copy(ψ)
        for v in vertices(pψ)
            pψ[v] = pψ[v] * onehot(only(s[v]) => bit_string[v])
        end

        bp_cache = build_bp_cache(copy(pψ))
        bp_cache = updatecache(bp_cache)

        f = scalar(Algorithm("bp"), bp_cache)
        p_bp = f*conj(f)

        # max_configuration_size = 12 corresponds to the sizes to consist of a heavy hexagon
        f = scalar(Algorithm("loopcorrections"), bp_cache; max_configuration_size=12)
        p_loop = f*conj(f)

        bmpsc = build_boundarymps_cache(copy(pψ), maxlinkdim(pψ)*2)
        bmpsc = updatecache(bmpsc)
        f = scalar(bmpsc)
        p_bmps = f*conj(f)

        println("p_bp = $p_bp")
        println("p_loop = $p_loop")
        println("p_bmps = $p_bmps")

        push!(probabilities_bp, p_bp)
        push!(probabilities_loop, p_loop)
        push!(probabilities_bmps, p_bmps)
    end
end

main()