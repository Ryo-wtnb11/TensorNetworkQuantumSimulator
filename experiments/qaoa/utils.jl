using NamedGraphs: NamedGraph
using NamedGraphs.GraphsExtensions: degree, neighbors

function parse_heavyhex_isingspinglass_setting_file(filename)
    txt = read(filename, String)
    txt = replace(txt, r"\btrue\b" => "true", r"\bfalse\b" => "false")
    txt = replace(txt, "{" => "Dict(", "}" => ")")
    txt = replace(txt, ":" => "=>")
    data = eval(Meta.parse(txt))
    dv_group, dij_group, dijk_group = data
    dv = merge(dv_group...)
    dij = merge(dij_group...)
    dijk = merge(dijk_group...)
    return dv, dij, dijk
end

function center_of_three_vertices(g::NamedGraph)
    center_of_three_neighbors_ = [v for v in vertices(g) if degree(g, v) == 2]
    center_of_three_neighbors = [v for v in center_of_three_neighbors_ if first(first(neighbors(g, v))) == first(last(neighbors(g, v))) || last(first(neighbors(g, v))) == last(last(neighbors(g, v)))]
    return center_of_three_neighbors
end

function heavyhex_isingspinglass_random_setting(g::NamedGraph)
    dv = Dict(v => rand([-1, 1]) for v in vertices(g))
    dij = Dict(e => rand([-1, 1]) for e in edges(g))
    center_of_three_neighbors = center_of_three_vertices(g)
    dijk = Dict(v => rand([-1, 1]) for v in center_of_three_neighbors)
    return dv, dij, dijk
end