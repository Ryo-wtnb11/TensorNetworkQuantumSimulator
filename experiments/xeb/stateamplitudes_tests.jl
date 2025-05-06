using Random
using Statistics
using Dictionaries

using TensorNetworkQuantumSimulator
const TN = TensorNetworkQuantumSimulator

using ITensorNetworks: ITensorNetworks, siteinds, Algorithm
const ITN = ITensorNetworks

using NamedGraphs: NamedEdge
using NamedGraphs.GraphsExtensions: add_vertices, add_edges

using Statistics

using ITensors

using StatsBase

Random.seed!(1734)

using Serialization

function main()

    g = TN.heavy_hexagonal_lattice(2,2)
    s = siteinds("S=1/2", g)
    ψ = ITN.random_tensornetwork(ComplexF64, s; link_space = 8)
    ψ, ψψ = normalize(ψ)

    n = 1 # You can change this to sample multiple amplitudes

    set_global_boundarymps_update_kwargs!(
        message_update_kwargs = (; niters = 25, tolerance = 1e-12),
    )

    bit_strings = BitVector[]
    probabilities = Float64[]

    seen = Set{BitVector}()

    for sample=1:n
        @show sample
        bit_string_info = [rand([1, 2]) for v in vertices(ψ)]
        bit_string = Dictionary(vertices(ψ), bit_string_info)

        bit = BitVector(bit_string_info .== 2)

        if bit in seen
            continue
        end
        push!(seen, bit)

        pψ = copy(ψ)
        for v in vertices(pψ)
            pψ[v] = pψ[v] * onehot(only(siteinds(pψ)[v]) => bit_string[v])
        end

        t = time()
        bmpsc = build_boundarymps_cache(copy(pψ), maxlinkdim(pψ))
        f = scalar(bmpsc)
        @info scalar(bmpsc)
        t = time() - t

        @which scalar(bmpsc)


        p = f * conj(f)
        push!(probabilities, real(p))
        push!(bit_strings, bit)
    end
end

main()