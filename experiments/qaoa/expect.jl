using NamedGraphs: src, dst, neighbors
include("utils.jl")

function expect_isingspinglass_without_zzz(g, ψ, dv, dij)
    σz = [real(expect(ψ, ("Z", [v], dv[v]); alg = "bp")) for v in vertices(g)]
    σzz = [real(expect(ψ, ("ZZ", [src(e), dst(e)], dij[e]); alg = "bp")) for e in edges(g)]
    σzzz = [real(expect(ψ, ("ZZZ", [first(neighbors(g, v)), v, last(neighbors(g, v))], dijk[v]); alg = "bp")) for v in center_of_three_neighbors]
    expect_val = sum(σz) + sum(σzz) + sum(σzzz)
    return expect_val
end

function expect_isingspinglass(g, ψ, dv, dij, dijk)
    center_of_three_neighbors = center_of_three_vertices(g)
    σz = [real(expect(ψ, ("Z", [v], dv[v]); alg = "bp")) for v in vertices(g)]
    σzz = [real(expect(ψ, ("ZZ", [src(e), dst(e)], dij[e]); alg = "bp")) for e in edges(g)]
    σzzz = [real(expect(ψ, ("ZZZ", [first(neighbors(g, v)), v, last(neighbors(g, v))], dijk[v]); alg = "bp")) for v in center_of_three_neighbors]
    expect_val = sum(σz) + sum(σzz) + sum(σzzz)
    return expect_val
end