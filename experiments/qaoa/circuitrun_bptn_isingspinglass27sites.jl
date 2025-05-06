using Downloads
using DelimitedFiles
using Random
using Serialization
using LinearAlgebra
Random.seed!(1234)

using TensorNetworkQuantumSimulator
const TN = TensorNetworkQuantumSimulator

using ITensorNetworks: ITensorNetworks, siteinds, Algorithm
const ITN = ITensorNetworks

using ITensors

using NamedGraphs: NamedEdge

using Statistics
using StatsBase

include("circuit_generators.jl")
include("utils.jl")
include("expect.jl")

# BLAS.set_num_threads(min(4, Sys.CPU_THREADS))
println("Julia is using "*string(Threads.nthreads()))
println("BLAS is using "*string(BLAS.get_num_threads()))
@show BLAS.get_config()

MAXDIM = 64

function main()
    # g = heavy_hexagonal_lattice_27sites()
    g = heavy_hexagonal_lattice_16sites()
    # initialize the random coefficients
    # local_file = Downloads.download("https://raw.githubusercontent.com/lanl/QAOA_vs_QA/main/problem_instances/ibm_geneva_0.txt")
    local_file = Downloads.download("https://raw.githubusercontent.com/lanl/QAOA_vs_QA/main/problem_instances/ibmq_guadalupe_0.txt")
    input_dv, input_dij, input_dijk = parse_heavyhex_isingspinglass_setting_file(local_file)
    dict_int_vert = deserialize("/Users/ryo/work/projects/QuantumAlgorithmTN/TensorNetworkQuantumSimulator/experiments/qaoa/input/index_namedvertex16sites.jld")

    center_of_three_neighbors = center_of_three_vertices(g)

    dv = Dict(dict_int_vert[key] => element for (key, element) in input_dv)
    dij = Dict((NamedEdge(dict_int_vert[first(key)] => dict_int_vert[last(key)]) in edges(g) ? reNamedEdge(dict_int_vert[first(key)] => dict_int_vert[last(key)]) : reverse(NamedEdge(dict_int_vert[first(key)] => dict_int_vert[last(key)]))) => element for (key, element) in input_dij)
    dijk = Dict(
        first(Tuple(filter(v -> v in center_of_three_neighbors, (dict_int_vert[i], dict_int_vert[j], dict_int_vert[k])))) => element
        for ((i, j, k), element) in input_dijk
    )

    # initialize the ITensorNetwork
    sites = siteinds("S=1/2", g)
    ψ = ITN.ITensorNetwork(v -> "↑", sites)
    ψψ = build_bp_cache(ψ)

    maxdim, cutoff = MAXDIM, 1e-14
    apply_kwargs = (; maxdim, cutoff, normalize = true)
    #Parameters for BP, as the graph is not a tree (it has loops), we need to specify these
    set_global_bp_update_kwargs!(;
        maxiter = 30,
        tol = 1e-10,
        message_update_kwargs = (;
            message_update_function = ms -> make_eigs_real.(ITN.default_message_update(ms))
        ),
    )

    no_trotter_steps = 5
    δt = 0.04

    hadamard_gates = [("H", [v]) for v in vertices(g)]
    ψ, ψψ, errors = apply(hadamard_gates, ψ, ψψ; apply_kwargs, verbose = false);

    for l = 1:no_trotter_steps
        layers = heavyhex_isingspinglass_circuit_layer(g, δt, δt, dv, dij, dijk)
        for layer in layers
            ψ, ψψ, errors = apply(layer, ψ, ψψ; apply_kwargs, verbose = false);
        end
    end

    expect_val = expect_isingspinglass(g, ψ, dv, dij, dijk)
    println("Expectation value is $expect_val")
end

main()