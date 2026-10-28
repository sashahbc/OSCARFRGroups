# spheres with a marking from a sphere group
#
# A marked sphere ("spider") is a point in Teichmueller space: a Delaunay
# triangulation of a set of points on P1, together with a marking
# homomorphism from a free group on a spanning-tree complement to a
# SphereGroup. Ported from gap/markedsphere.gd + .gi (infrastructure layer
# only -- see TODO for the deferred branched-covering machinery).

export MarkedSphere
export marking
export new_marked_sphere
export spanning_tree_boundary
export sphere_model
export tree_cost
export tree_edges
export triangulation
export vertices_of_marked_sphere
export wiggled_marked_sphere

struct MarkedSphere
  model::SphereGroup
  cut::SphereTriangulation
  group::GapObj
  intree::Vector{Bool}
  treeedge::Vector{TEdge}
  treecost::Float64
  marking::GapObj
end

triangulation(s::MarkedSphere) = s.cut
Oscar.group(s::MarkedSphere) = s.group
sphere_model(s::MarkedSphere) = s.model
marking(s::MarkedSphere) = s.marking
tree_cost(s::MarkedSphere) = s.treecost
tree_edges(s::MarkedSphere) = s.treeedge
in_tree(s::MarkedSphere) = s.intree
in_tree(s::MarkedSphere, e::TEdge) = in_tree(s)[index(e)]

Base.show(io::IO, s::MarkedSphere) = print(
  io, "<marked sphere on ", triangulation(s), " marked by ",
  String(GAP.Globals.String(marking(s))), ">",
)

"""the positions of the (non-fake) vertices a marked sphere lies on."""
vertices_of_marked_sphere(s::MarkedSphere) =
  [pos(v) for v in vertices(triangulation(s)) if !is_fake(v)]

################################################################
# spanning tree (triangulations.gi has no analogue -- generic Prim's
# algorithm on the triangulation's vertex/edge graph, weighted by
# edge_length; this is ordinary graph theory, not IMG-specific, so it does
# not transliterate GAP's MINSPANTREE@, just matches its contract: one
# directed TEdge per tree edge (oriented from the already-grown part of the
# tree to the newly attached vertex), plus the tree's total length.)
function _minimum_spanning_tree(cut::SphereTriangulation)
  nv = length(vertices(cut))
  visited = falses(nv)
  visited[1] = true
  dist = fill(Inf, nv)
  bestedge = Vector{Union{Nothing,TEdge}}(nothing, nv)
  for e in neighbors(vertices(cut)[1])
    j = index(edge_to(e))
    c = edge_length(e)
    c < dist[j] && (dist[j] = c; bestedge[j] = e)
  end

  tree = TEdge[]
  cost = 0.0
  for _ in 2:nv
    j = 0
    best = Inf
    for i in 1:nv
      if !visited[i] && dist[i] < best
        best = dist[i]
        j = i
      end
    end
    j == 0 && error("triangulation's vertex graph is disconnected")
    push!(tree, bestedge[j])
    cost += dist[j]
    visited[j] = true
    for e in neighbors(vertices(cut)[j])
      k = index(edge_to(e))
      if !visited[k] && edge_length(e) < dist[k]
        dist[k] = edge_length(e)
        bestedge[k] = e
      end
    end
  end
  tree, cost
end

"""edges traversed when walking once around the spanning tree with it on
our right, starting at the vertex with the highest index."""
function spanning_tree_boundary(cut::SphereTriangulation, intree::AbstractVector{Bool})
  n = length(vertices(cut))
  idx = findfirst(e -> intree[index(e)] && index(edge_from(e)) == n, edges(cut))
  idx === nothing && error("no tree edge found starting at vertex $n")
  e = edges(cut)[idx]
  first_edge = e
  result = TEdge[]
  while true
    push!(result, e)
    e = opposite(e)
    while true
      e = nextedge(opposite(e))
      intree[index(e)] && break
    end
    e === first_edge && break
  end
  result
end
spanning_tree_boundary(s::MarkedSphere) =
  spanning_tree_boundary(triangulation(s), in_tree(s))

################################################################
# construction (markedsphere.gi:111-206)

"""
    new_marked_sphere(points::AbstractVector{<:P1Point}, model=sphere_group(length(points)); quality=100.0)

The marked sphere ("spider") with identity marking on `points`: a Delaunay
triangulation of `points`, with a spanning tree of its vertex graph marked
by a free group on the tree complement, mapped into `model` via a marking
homomorphism."""
function new_marked_sphere(
  points::Vector{P1Point{T}}, model::SphereGroup=sphere_group(length(points));
  quality::Real=100.0,
) where {T}
  n = length(points)
  sphere_rank(model) == n ||
    error("model must have the same rank as the number of points")

  cut = delaunay_triangulation(points; quality)
  tree, treecost = _minimum_spanning_tree(cut)
  sort!(tree; by=edge_length)

  # start by a free group on the tree edges; by convention, each tree edge's
  # positively-oriented generator goes in the direction Prim's algorithm
  # discovered it.
  g0 = GAP.Globals.FreeGroup(length(tree))
  g0gens = GAP.Globals.GeneratorsOfGroup(g0)

  for e in edges(cut)
    group_element!(e, one(g0))
  end
  intree = falses(length(edges(cut)))
  for i in 1:length(tree)
    x = g0gens[i]
    group_element!(tree[i], x)
    group_element!(opposite(tree[i]), inv(x))
    intree[index(tree[i])] = true
    intree[index(opposite(tree[i]))] = true
  end

  p = GAP.Globals.PresentationFpGroup(g0, 0)
  GAP.Globals.TzOptions(p).protected = length(tree)
  GAP.Globals.TzInitGeneratorImages(p)

  # add relators saying the cycle around each fake vertex is trivial
  for v in vertices(cut)
    if is_fake(v)
      relator = GAP.Globals.Product(
        GapObj(reverse([group_element(e) for e in neighbors(v)]))
      )
      GAP.Globals.AddRelator(p, relator)
    end
  end

  # eliminate every tree-edge generator that a fake-vertex relator allows;
  # exactly length(tree)-(n-1) of them (one per fake vertex) can go.
  for pg in GAP.Globals.GeneratorsOfPresentation(p)
    GAP.Globals.TzEliminate(p, pg)
  end

  treegroup = GAP.Globals.FreeGroup(n - 1)
  treegroupgens = GAP.Globals.GeneratorsOfGroup(treegroup)

  for e in edges(cut)
    group_element!(e, one(treegroup))
  end
  images = GAP.Globals.TzImagesOldGens(p)
  pgens = GAP.Globals.GeneratorsOfPresentation(p)
  for i in 1:length(tree)
    x = GAP.Globals.MappedWord(images[i], pgens, treegroupgens)
    group_element!(tree[i], x)
    group_element!(opposite(tree[i]), inv(x))
  end

  treeedge = [
    tree[Int(GAP.Globals.TietzeWordAbstractWord(w)[1])] for
    w in GAP.Globals.TzPreImagesNewGens(p)
  ]

  # read off, for each real vertex, the loop word around it, to build the
  # marking (treegroup -> model)
  ordering = Int[]
  source = Dict{Int,GapObj}()
  for e in spanning_tree_boundary(cut, intree)
    if !is_fake(edge_from(e))
      i = index(edge_from(e))
      if !haskey(source, i)
        source[i] = one(treegroup)
        push!(ordering, i)
      end
      source[i] = GAP.Globals.LeftQuotient(group_element(e), source[i])
    end
  end

  fmodel = GAP.Globals.FreeGroupOfFpGroup(GapObj(model))
  gmodel = GAP.Globals.GeneratorsOfGroup(fmodel)
  conjugated(ord, i) =
    gmodel[ord[i]]^GAP.Globals.Inverse(
      GAP.Globals.Product(
        GapObj([gmodel[j] for j in sort(collect(setdiff(1:ord[i], ord[1:i])))])
      ),
    )

  image = GAP.Globals.GroupHomomorphismByImages(
    fmodel, fmodel, GapObj([gmodel[i] for i in ordering]),
    GapObj([conjugated(ordering, i) for i in 1:n]),
  )
  modelordering = sphere_ordering(model)
  image2 = GAP.Globals.GroupHomomorphismByImages(
    fmodel, fmodel, GapObj([conjugated(modelordering, i) for i in 1:n]),
    GapObj([gmodel[i] for i in modelordering]),
  )
  image = image * image2

  fam = GAP.Globals.FamilyObj(one(GapObj(model)))
  marking = GAP.Globals.GroupHomomorphismByImages(
    treegroup, GapObj(model), GapObj([source[i] for i in 1:n]),
    GapObj([
      GAP.Globals.ElementOfSphereGroup(fam, GAP.Globals.Image(image, gmodel[x])) for
      x in 1:n
    ]),
  )

  MarkedSphere(model, cut, treegroup, intree, treeedge, treecost, marking)
end

################################################################
# wiggling (markedsphere.gi:215-235)

"""rebuild the marked sphere s with its underlying triangulation
wiggled (see [`wiggled_triangulation`](@ref)); the marking is unaffected."""
function wiggled_marked_sphere(s::MarkedSphere, movement)
  newcut = wiggled_triangulation(triangulation(s), movement)
  newtreeedge = [edges(newcut)[index(e)] for e in tree_edges(s)]
  MarkedSphere(
    sphere_model(s), newcut, Oscar.group(s), in_tree(s), newtreeedge, tree_cost(s),
    marking(s),
  )
end

shallow_copy(s::MarkedSphere) = wiggled_marked_sphere(s, nothing)
