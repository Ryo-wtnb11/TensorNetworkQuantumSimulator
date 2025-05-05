function SimpleGraphAlgorithms.edge_color(g::AbstractGraph, k::Int64)
    pg, vs = position_graph(g), collect(vertices(g))
    ec_dict = edge_color(UG(pg), k)
    # returns k vectors which contain the colored/commuting edges
    return [
        [(vs[first(first(e))], vs[last(first(e))]) for e in ec_dict if last(e) == i] for
        i = 1:k
    ]
end

"""Create heavy-hex lattice geometry"""
function heavy_hexagonal_lattice(nx::Int64, ny::Int64)
    g = named_hexagonal_lattice_graph(nx, ny)
    # create some space for inserting the new vertices
    g = rename_vertices(v -> (2 * first(v) - 1, 2 * last(v) - 1), g)
    for e in edges(g)
        vsrc, vdst = src(e), dst(e)
        v_new = ((first(vsrc) + first(vdst)) / 2, (last(vsrc) + last(vdst)) / 2)
        g = add_vertex(g, v_new)
        g = rem_edge(g, e)
        g = add_edges(g, [NamedEdge(vsrc => v_new), NamedEdge(v_new => vdst)])
    end
    return g
end

function lieb_lattice(nx::Int64, ny::Int64; periodic = false)
    @assert (!periodic && isodd(nx) && isodd(ny)) || (periodic && iseven(nx) && iseven(ny))
    g = named_grid((nx, ny); periodic)
    for v in vertices(g)
        if iseven(first(v)) && iseven(last(v))
            g = rem_vertex(g, v)
        end
        if iseven(first(v)) && iseven(last(v))
            g = rem_vertex(g, v)
        end
    end
    return g

end



function topologytograph(topology)
    # TODO: adapt this to named graphs with non-integer labels
    # find number of vertices
    nq = maximum(maximum.(topology))
    adjm = zeros(Int, nq, nq)
    for (ii, jj) in topology
        adjm[ii, jj] = adjm[jj, ii] = 1
    end
    return NamedGraph(SimpleGraph(adjm))
end


function graphtotopology(g)
    return [[edge.src, edge.dst] for edge in edges(g)]
end

function NamedGraphs.GraphsExtensions.rem_vertex(bpc::AbstractBeliefPropagationCache, v)
    return rem_vertices(bpc, [v])
end

function NamedGraphs.GraphsExtensions.rem_vertices(bpc::BeliefPropagationCache, vs::Vector)
    pg = partitioned_tensornetwork(bpc)
    pg = rem_vertices(pg, vs)
    return BeliefPropagationCache(pg, messages(bpc))
end

function NamedGraphs.GraphsExtensions.rem_vertices(bmpsc::BoundaryMPSCache, vs::Vector)
    bpc = bp_cache(bmpsc)
    bpc = rem_vertices(bpc, vs)
    return BoundaryMPSCache(bpc, ppg(bmpsc), maximum_virtual_dimension(bmpsc))
end

# for sycamore architecture
function rem_triangular_block(g::NamedGraph, nxny::Tuple, p1::Tuple, p2::Tuple, delete_side::Symbol)
    nx, ny = nxny
    g = copy(g)
    x1, y1 = p1
    x2, y2 = p2
    dx = x2 - x1
    dy = y2 - y1
    for y in 1:ny, x in 1:nx
        determinant = dy * (x - x1) - dx * (y - y1)
        remove_vertex = (delete_side == :left && determinant > 0) ||
                        (delete_side == :right && determinant < 0) ||
                        (delete_side == :on_line && determinant == 0)
        if remove_vertex && has_vertex(g, (x, y))
            g = rem_vertex(g, (x, y))
        end
    end
    return g
end

function rem_triangular_blocks(g::NamedGraph, nxny::Tuple, blocks::Vector{Tuple{Tuple{Int, Int}, Tuple{Int, Int}, Symbol}})
    nx, ny = nxny
    g = copy(g)
    nx, ny = last(vertices(g))
    for (p1, p2, delete_side) in blocks
        g = rem_triangular_block(g, (nx, ny), p1, p2, delete_side)
    end
    return g
end

# Fig .1a of https://www.nature.com/articles/s41586-019-1666-5
function sycamore_53_qubit_grid(;patch::Int64=0)
    nx, ny = 10, 10
    g = named_grid((nx, ny))
    blocks = [
        ((1, 6), (6, 1), :left),
        ((7, 1), (10, 4), :left),
        ((1, 6), (5, 10), :right),
        ((5, 10), (10, 5), :right),
    ]
    g = rem_triangular_blocks(g, (nx, ny), blocks)
    g = rem_vertex(g, (4, 3))
    if patch == 1
        g = rem_triangular_block(g, (nx, ny), (4, 4), (7, 7), :left)
    elseif patch == 2
        g = rem_triangular_block(g, (nx, ny), (5, 4), (8, 7), :right)
    end
    return g
end

# Fig .4b of https://www.nature.com/articles/s41586-024-07998-6
function sycamore_66_qubit_grid(;reduced = false)
    nx, ny = 11, 11
    g = named_grid((nx, ny))
    blocks = [
        ((1, 5), (5, 1), :left),
        ((6, 1), (11, 6), :left),
        ((1, 6), (6, 11), :right),
        ((6, 11), (11, 6), :right),
    ]
    g = rem_triangular_blocks(g, (nx, ny), blocks)
    if reduced
        g = rem_triangular_block(g, (nx, ny), (1, 6), (6, 1), :left)
        g = rem_vertex(g, [(1, 6), (4, 3),(6, 1), (6, 11), (11, 6)])
    end
    return g
end

# Fig .4a of https://www.nature.com/articles/s41586-024-07998-6
function sycamore_67_qubit_grid()
    nx, ny = 11, 11
    g = named_grid((nx, ny))
    blocks = [
        ((1, 6), (6, 1), :left)
        ((7, 1), (11, 5), :left)
        ((1, 6), (6, 11), :right)
        ((7, 11), (11, 7), :right)
    ]
    g = rem_triangular_blocks(g, (nx, ny), blocks)
    g = rem_vertex(g, [(7, 1), (7, 2), (10, 7), (11, 7)])
    return g
end

# Fig. 4-6 of https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-024-07998-6/MediaObjects/41586_2024_7998_MOESM1_ESM.pdf
function sycamore_70_qubit_grid()
    nx, ny = 11, 11
    g = named_grid((nx, ny))
    blocks = [
        (g, (1, 5), (5, 1), :left)
        (g, (7, 1), (11, 5), :left)
        (g, (1, 6), (6, 11), :right)
        (g, (6, 11), (11, 6), :right)
    ]
    g = rem_vertex(g, (6, 1))
    return g
end
