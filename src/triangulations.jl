# sphere triangulations
#
# A half-edge (DCEL) representation of a triangulated sphere: vertices,
# edges (in opposite pairs, each with a face on its left), and faces
# (always triangles). Ported from gap/triangulations.gd + .gi.
#
# Unlike GAP, derived quantities (edge midpoint/length/map, face
# centre/radius/barycentre) are NOT cached -- they're recomputed on the fly
# from vertex positions every time. This drops GAP's HasXXX/SetXXX/
# ResetFilterObj bookkeeping entirely; the only genuinely stored fields are
# the topology pointers (from/left/nextedge/prevopp/opposite, neighbour)
# plus each vertex's position/is_fake flag and each edge's optional
# GroupElement.

export TVertex, TEdge, TFace, SphereTriangulation
export neighbour, neighbours, valency, edge_from, edge_to, edge_left, edge_right
export nextedge, prevopp, opposite
export edge_map, edge_length, from_pos, to_pos, pos, radius, centre
export group_element, group_element!, has_group_element, is_fake
export edge_path, locate, locate_face_in_triangulation, locate_in_triangulation
export InFace, OnEdge, AtVertex, face_of, located_object
export flip_edge!, check_triangulation, fix_delaunay!
export add_to_triangulation!, remove_from_triangulation!
export delaunay_triangulation, equidistributed_p1_points
export wiggled_triangulation, shallow_copy
export closest_face, closest_faces, closest_vertex, closest_vertices

################################################################
# types

abstract type TriangulationObject end

mutable struct TVertex <: TriangulationObject
    index::Int
    neighbour::Union{Nothing,TriangulationObject}   # some edge starting at this vertex
    pos::Union{Nothing,P1Point}
    is_fake::Bool
    TVertex(index::Integer) = new(index, nothing, nothing, false)
end

mutable struct TEdge <: TriangulationObject
    index::Int
    from::Union{Nothing,TriangulationObject}        # TVertex
    left::Union{Nothing,TriangulationObject}         # TFace
    nextedge::Union{Nothing,TriangulationObject}     # TEdge: next edge around left face
    prevopp::Union{Nothing,TriangulationObject}      # TEdge: opposite(previous edge around left face)
    opposite::Union{Nothing,TriangulationObject}     # TEdge: reversed edge (fixed pairing, set once)
    groupelement::Any                                 # nothing, or set explicitly
    TEdge(index::Integer) = new(index, nothing, nothing, nothing, nothing, nothing, nothing)
end

mutable struct TFace <: TriangulationObject
    index::Int
    neighbour::Union{Nothing,TriangulationObject}   # some edge with this face on its left
    TFace(index::Integer) = new(index, nothing)
end

struct SphereTriangulation
    v::Vector{TVertex}
    e::Vector{TEdge}
    f::Vector{TFace}
end

################################################################
# navigation & display (triangulations.gi:15-120)

Base.show(io::IO, v::TVertex) = print(io, "<vertex ", v.index, " ", [e.index for e in neighbours(v)], ">")
Base.show(io::IO, e::TEdge) = print(io, "<edge ", e.index, " ", [edge_from(e).index, edge_to(e).index], ">")
Base.show(io::IO, f::TFace) = print(io, "<face ", f.index, " ", [e.index for e in neighbours(f)], ">")
Base.show(io::IO, t::SphereTriangulation) =
    print(io, "<triangulation with ", length(t.v), " vertices, ", length(t.e), " edges and ", length(t.f), " faces>")

edge_from(e::TEdge) = e.from
edge_left(e::TEdge) = e.left
nextedge(e::TEdge) = e.nextedge
prevopp(e::TEdge) = e.prevopp
opposite(e::TEdge) = e.opposite
neighbour(v::TVertex) = v.neighbour
neighbour(f::TFace) = f.neighbour
is_fake(v::TVertex) = v.is_fake

edge_to(e::TEdge) = edge_from(opposite(e))
edge_right(e::TEdge) = edge_left(opposite(e))

has_group_element(e::TEdge) = e.groupelement !== nothing
function group_element(e::TEdge)
    e.groupelement !== nothing && return e.groupelement
    oe = opposite(e)
    oe.groupelement !== nothing && return inv(oe.groupelement)
    error("no GroupElement set for this edge or its opposite")
end
group_element!(e::TEdge, g) = (e.groupelement = g)

function neighbours(v::TVertex, e0::TEdge = neighbour(v))
    n = TEdge[]
    e = e0
    while true
        push!(n, e)
        e = prevopp(e)
        e == e0 && break
    end
    n
end
valency(v::TVertex) = length(neighbours(v))

function neighbours(f::TFace, e0::TEdge = neighbour(f))
    n = TEdge[]
    e = e0
    while true
        push!(n, e)
        e = nextedge(e)
        e == e0 && break
    end
    n
end
function valency(f::TFace)
    n = neighbours(f)
    length(n) == 3 || error("found a face with $(length(n)) != 3 adjacent edges")
    3
end

pos(v::TVertex) = v.pos
from_pos(e::TEdge) = pos(edge_from(e))
to_pos(e::TEdge) = pos(edge_to(e))
pos(e::TEdge) = barycentre(from_pos(e), to_pos(e))
edge_map(e::TEdge) = moebius_path(from_pos(e), to_pos(e))
edge_length(e::TEdge) = distance(from_pos(e), to_pos(e))

"""the circumcentre and circumradius of face f's 3 vertices (as returned
jointly by p1points.jl's `circumcircle`, to avoid computing it twice)."""
circumcircle_data(f::TFace) = circumcircle(from_pos.(neighbours(f))...)
centre(f::TFace) = circumcircle_data(f)[1]
radius(f::TFace) = circumcircle_data(f)[2]
pos(f::TFace) = barycentre(from_pos.(neighbours(f)))

################################################################
# edge_path (triangulations.gi:175-198)

function edge_path(t::SphereTriangulation, f0::TFace, f1::TFace)
    stuck = length(t.f)
    path = TEdge[]
    while f0 != f1
        stuck -= 1
        stuck < 0 && error("stuck in a loop trying to find a path")
        n = neighbours(f0)
        ymin = ycoord(preimage(edge_map(n[1]), pos(f1)))
        emin = n[1]
        for e in n[2:end]
            y = ycoord(preimage(edge_map(e), pos(f1)))
            if y < ymin
                ymin = y
                emin = e
            end
        end
        push!(path, emin)
        f0 = edge_right(emin)
    end
    path
end

################################################################
# point location (triangulations.gi:200-291)

abstract type LocateResult end
struct InFace <: LocateResult
    face::TFace
    yc::Vector
end
struct OnEdge <: LocateResult
    face::TFace
    edge::TEdge
    t::Any
end
struct AtVertex <: LocateResult
    face::TFace
    edge_in::TEdge
    edge_out::TEdge
    vertex::TVertex
end

face_of(r::LocateResult) = r.face
located_object(r::InFace) = r.face
located_object(r::OnEdge) = r.edge
located_object(r::AtVertex) = r.vertex

"""walk from face f0 towards p, returning an InFace/OnEdge/AtVertex result.

Unlike GAP (which uses numeric epsilon thresholds appropriate to its default
IEEE754 backend), this uses exact/rigorous sign tests (`>0`/`<0`) suited to
the exact (QQBarField) and ball-arithmetic (AcbField) backends this port
targets: a y-ratio that's neither provably positive nor provably negative is
treated as the boundary case."""
function locate(t::SphereTriangulation, f0::TFace, p::P1Point)
    stuck = length(t.f)
    while true
        stuck -= 1
        stuck < 0 && error("stuck in a loop trying to locate a face")

        n = neighbours(f0)
        baryc = [preimage(edge_map(e), p) for e in n]
        yc = ycoord.(baryc)
        pos_ = [y > 0 for y in yc]

        all(pos_) && return InFace(f0, yc)

        neg = findfirst(i -> yc[i] < 0, 1:3)
        if neg !== nothing
            f0 = edge_right(n[neg])
            continue
        end

        zeros_ = findall(i -> !pos_[i], 1:3)
        if length(zeros_) == 1
            i = zeros_[1]
            return OnEdge(f0, n[i], real(zcoord(baryc[i])))
        elseif length(zeros_) == 2
            i, j = zeros_
            (1 + i % 3) == j || ((i, j) = (j, i))
            return AtVertex(f0, n[i], n[j], edge_to(n[i]))
        else
            error("degenerate triangle (flat angle) while locating a point")
        end
    end
end

locate_face_in_triangulation(t::SphereTriangulation, p::P1Point) = face_of(locate(t, t.f[1], p))
locate_face_in_triangulation(t::SphereTriangulation, seed::TVertex, p::P1Point) = face_of(locate(t, edge_left(neighbour(seed)), p))
locate_face_in_triangulation(t::SphereTriangulation, seed::TEdge, p::P1Point) = face_of(locate(t, edge_left(seed), p))
locate_face_in_triangulation(t::SphereTriangulation, seed::TFace, p::P1Point) = face_of(locate(t, seed, p))

locate_in_triangulation(t::SphereTriangulation, p::P1Point) = located_object(locate(t, t.f[1], p))
locate_in_triangulation(t::SphereTriangulation, seed::TVertex, p::P1Point) = located_object(locate(t, edge_left(neighbour(seed)), p))
locate_in_triangulation(t::SphereTriangulation, seed::TEdge, p::P1Point) = located_object(locate(t, edge_left(seed), p))
locate_in_triangulation(t::SphereTriangulation, seed::TFace, p::P1Point) = located_object(locate(t, seed, p))

################################################################
# edge flip, consistency check, Delaunay fixup (triangulations.gi:315-401)

y_ratio(a::P1Point,b::P1Point,c::P1Point,d::P1Point) = ycoord(x_ratio(a,b,c,d))

"""flip the edge e. If multi, recursively keep flipping (propagating from
the two new diagonals) until the neighbourhood is Delaunay; if not multi,
perform exactly one flip."""
function flip_edge!(e::TEdge, multi::Bool)
    f = opposite(e)
    a = edge_from(e); b = edge_from(f)
    bp = nextedge(e); pa = nextedge(bp)
    aq = nextedge(f); qb = nextedge(aq)
    p = edge_from(pa); q = edge_from(qb)

    (!multi || y_ratio(pos(p), pos(q), pos(a), pos(b)) > 0) || return

    paq = edge_left(f)
    qbp = edge_left(e)
    opa = opposite(pa); oqb = opposite(qb)

    e.from = p; e.nextedge = qb; e.prevopp = opposite(bp)
    f.from = q; f.nextedge = pa; f.prevopp = opposite(aq)

    qbp.neighbour = e
    paq.neighbour = f
    a.neighbour = aq
    b.neighbour = bp

    pa.prevopp = e; pa.nextedge = aq; pa.left = paq
    qb.prevopp = f; qb.nextedge = bp; qb.left = qbp
    aq.prevopp = opa; aq.nextedge = f
    bp.prevopp = oqb; bp.nextedge = e

    if has_group_element(e)
        pa.groupelement = inv(group_element(e)) * group_element(pa)
        opposite(pa).groupelement = inv(pa.groupelement)
        qb.groupelement = group_element(e) * group_element(qb)
        opposite(qb).groupelement = inv(qb.groupelement)
    end

    if multi
        flip_edge!(aq, true)
        flip_edge!(qb, true)
    end
    nothing
end

"""verify the DCEL invariants of t (and, unless nodelaunay, the Delaunay
condition too). Returns `true`, or a (message, offending objects) pair."""
function check_triangulation(t::SphereTriangulation; nodelaunay::Bool=false)
    bad = filter(v -> edge_from(neighbour(v)) !== v, t.v)
    isempty(bad) || return ("From(Neighbour(v)) != v: ", bad)

    bad = filter(e -> opposite(opposite(e)) !== e, t.e)
    isempty(bad) || return ("Opposite(Opposite(e)) != e: ", bad)

    bad = filter(e -> opposite(e) !== prevopp(nextedge(e)), t.e)
    isempty(bad) || return ("Opposite(e) != Prevopp(NextEdge(e)): ", bad)

    bad = filter(e -> e !== nextedge(nextedge(nextedge(e))), t.e)
    isempty(bad) || return ("NextEdge(NextEdge(NextEdge(e))) != e: ", bad)

    bad = filter(e -> edge_from(e) !== edge_from(prevopp(e)), t.e)
    isempty(bad) || return ("From(e) != From(Prevopp(e)): ", bad)

    bad = filter(e -> edge_left(e) !== edge_left(nextedge(e)), t.e)
    isempty(bad) || return ("Left(e) != Left(NextEdge(e)): ", bad)

    # Left(e) == Right(Opposite(e)) and From(e) == To(Opposite(e)) hold
    # tautologically here, since edge_right/edge_to are defined directly in
    # terms of `opposite` -- nothing to check.

    bad = filter(f -> edge_left(neighbour(f)) !== f, t.f)
    isempty(bad) || return ("Left(Neighbour(f)) != f: ", bad)

    if !nodelaunay
        bad = filter(t.e) do e
            y_ratio(to_pos(nextedge(e)), to_pos(nextedge(opposite(e))), from_pos(e), to_pos(e)) > 0
        end
        isempty(bad) || return ("Delaunay condition fails: ", bad)

        bad = filter(t.f) do f
            locate_in_triangulation(t, f, pos(f)) !== f
        end
        isempty(bad) || return ("Locate(f,Pos(f)) != f: ", bad)
    end

    true
end

function fix_delaunay!(t::SphereTriangulation)
    fixes = 0
    idle = false
    while !idle
        idle = true
        for e in t.e
            if y_ratio(to_pos(nextedge(e)), to_pos(nextedge(opposite(e))), from_pos(e), to_pos(e)) > 0
                flip_edge!(e, true)
                fixes += 1
                idle = false
            end
        end
    end
    fixes
end

################################################################
# add/remove a vertex (triangulations.gi:403-520)

"""add the P1 point p, located in face f0, to t. If delaunay, flip diagonals
around the new vertex as needed to preserve the Delaunay condition."""
function add_to_triangulation!(t::SphereTriangulation, f0::TFace, p::P1Point; delaunay::Bool=true)
    v = TVertex(length(t.v)+1)
    push!(t.v, v)
    v.pos = p

    oldn = TEdge[]
    local e2::TEdge
    for e0 in neighbours(f0)
        push!(oldn, e0)
        f = TFace(length(t.f)+1); push!(t.f, f)
        e1 = TEdge(length(t.e)+1); push!(t.e, e1)
        e2 = TEdge(length(t.e)+1); push!(t.e, e2)

        f.neighbour = e0
        e1.from = edge_to(e0)
        e2.from = v
        e1.left = f
        e2.left = f
        e0.left = f
        # Opposite(e0)'s Right becomes f automatically: edge_right is derived.

        e1.prevopp = opposite(e0)
        old_next = nextedge(e0)
        old_next.prevopp = e1
        e1.nextedge = e2
        e2.nextedge = e0
        e0.nextedge = e1

        if has_group_element(e0)
            ge = one(group_element(e0))
            e1.groupelement = ge
            e2.groupelement = ge
        end
    end
    v.neighbour = e2

    for e0 in oldn
        x = nextedge(nextedge(e0))
        y = prevopp(e0)
        x.opposite = y
        y.opposite = x
    end
    for e0 in oldn
        nextedge(nextedge(e0)).prevopp = opposite(nextedge(e0))
    end

    # recycle f0's slot with the just-created last face
    last_f = pop!(t.f)
    t.f[f0.index] = last_f
    last_f.index = f0.index

    if delaunay
        for e0 in oldn
            flip_edge!(e0, true)
        end
    end
    v
end
add_to_triangulation!(t::SphereTriangulation, p::P1Point; seed::Union{Nothing,TFace}=nothing, delaunay::Bool=true) =
    add_to_triangulation!(t, seed === nothing ? face_of(locate(t, t.f[1], p)) : seed, p; delaunay)

"""remove vertex v from t: flip edges as needed until v becomes trivalent,
then delete it and its 3 incident edges/2 faces."""
function remove_from_triangulation!(t::SphereTriangulation, v::TVertex)
    while valency(v) > 3
        e0 = neighbour(v)
        e1 = prevopp(e0)
        e2 = prevopp(e1)
        miny = Inf
        echosen = e1
        start = e0
        while true
            c, r = circumcircle(to_pos(e0), to_pos(e1), to_pos(e2))
            y = r^2 - distance(c, pos(v))^2
            if y < miny
                miny = y
                echosen = e1
            end
            e0, e1, e2 = e1, e2, prevopp(e2)
            e0 == start && break
        end
        flip_edge!(echosen, false)
    end

    e0 = neighbour(v); e1 = prevopp(e0); e2 = prevopp(e1)

    # deallocate faces Left(e1) and Left(e2)
    for e in (e1, e2)
        j = edge_left(e).index
        f = pop!(t.f)
        if j <= length(t.f)
            t.f[j] = f
            f.index = j
        end
    end

    # deallocate the 6 edges in/out of v
    for e in (e0, e1, e2, opposite(e0), opposite(e1), opposite(e2))
        j = e.index
        e_last = pop!(t.e)
        if j <= length(t.e)
            t.e[j] = e_last
            e_last.index = j
        end
    end

    # deallocate vertex v
    j = v.index
    v_last = pop!(t.v)
    if j <= length(t.v)
        t.v[j] = v_last
        v_last.index = j
    end

    e0 = nextedge(e0); e1 = nextedge(e1); e2 = nextedge(e2)

    f = edge_left(e0)
    f.neighbour = e0
    edge_from(e0).neighbour = e0; e0.left = f; e0.nextedge = e1; e0.prevopp = opposite(e2)
    edge_from(e1).neighbour = e1; e1.left = f; e1.nextedge = e2; e1.prevopp = opposite(e0)
    edge_from(e2).neighbour = e2; e2.left = f; e2.nextedge = e0; e2.prevopp = opposite(e1)
    nothing
end

################################################################
# Delaunay triangulation construction (triangulations.gi:522-624)

const _OCTAHEDRON_ADJACENCY = (
    (1,2,1,3),(2,6,4,13),(3,6,1,23),(4,5,7,17),(5,5,8,19),(6,1,7,4),
    (7,1,5,9),(8,3,6,18),(9,3,5,20),(10,4,3,14),(11,4,2,24),(12,2,3,10),
    (13,2,4,15),(14,3,3,12),(15,3,4,2),(16,6,6,8),(17,6,7,6),(18,1,6,16),
    (19,1,8,21),(20,4,5,7),(21,4,8,5),(22,5,2,11),(23,5,1,1),(24,2,2,22),
)

"""construct a Delaunay triangulation of the sphere with the given points as
vertices. If some points are aligned on a great circle, or all lie in a
hemisphere, extra ("fake", `is_fake(v)==true`) vertices are added to make an
initial simplicial octahedron possible. If quality is given, faces are
further refined (adding circumcentres, also marked fake) until every face's
circumradius/(shortest edge) ratio is at most quality."""
function delaunay_triangulation(points::Vector{P1Point{T}}; quality::Real=Inf) where T
    points = copy(points)
    n = length(points)
    if n == 0
        points = [P1Inf(T)]
        n = 1
    end

    F = field(points[n])
    dist_n = [distance(points[n], x) for x in points]
    order = [n]
    dev90 = [abs(v - pi/2) for v in dist_n]
    ibest = argmin(dev90)

    if dev90[ibest] >= pi/6 # all points roughly aligned with points[n]
        ifar = argmax(dist_n)
        if dist_n[ifar] < pi/2 # actually all points are close to points[n]
            push!(points, antipode(points[n]))
            push!(order, length(points))
        else
            push!(order, ifar)
        end
        mp = moebius_map(points[order[2]], points[n]) # 0 ↦ points[n], ∞ ↦ points[order[2]]
        for p in (P1one(F), P1Point(F(im)), P1Point(-one(F)), P1Point(F(-im)))
            push!(points, image(mp, p))
            push!(order, length(points))
        end
    else # points[ibest] is roughly at 90 degrees from points[n]
        mp = moebius_map(antipode(points[n]), points[n], points[ibest]) # 0↦points[n], 1↦points[ibest], ∞↦antipode(points[n])
        invmp = inv(mp)
        mapped = [image(invmp, x) for x in points]
        for p in (P1Inf(F), P1one(F), P1Point(F(im)), P1Point(-one(F)), P1Point(F(-im)))
            dd = [distance(x, p) for x in mapped]
            j = argmin(dd)
            if dd[j] >= pi/6
                push!(points, image(mp, p))
                push!(order, length(points))
            else
                push!(order, j)
            end
        end
    end
    allunique(order) || error("delaunay_triangulation could not create initial octahedron")
    append!(order, setdiff(1:n, order))

    # build the initial octahedron
    v = [TVertex(order[i]) for i in 1:6]
    e = [TEdge(i) for i in 1:24]
    f = [TFace(i) for i in 1:8]
    for i in 1:6
        v[i].pos = points[order[i]]
    end
    for i in 1:24
        e[i].opposite = e[iseven(i) ? i-1 : i+1]
    end
    for (ei,vi,fi,ni) in _OCTAHEDRON_ADJACENCY
        e[ei].from = v[vi]
        v[vi].neighbour = e[ei]
        e[ei].left = f[fi]
        f[fi].neighbour = e[ei]
        e[ei].nextedge = e[ni]
        e[ni].prevopp = e[ei].opposite
    end
    t = SphereTriangulation(v, e, f)

    # insert the remaining points
    for i in 7:length(order)
        p = points[order[i]]
        result = locate(t, t.f[1], p)
        result isa AtVertex && error("Two vertices coincide: $p and $(result.vertex)")
        add_to_triangulation!(t, face_of(result), p; delaunay=true)
        t.v[i].index = order[i]
    end

    # restore the original point ordering
    vcopy = copy(t.v)
    for k in 1:length(order)
        t.v[order[k]] = vcopy[k]
    end

    # quality refinement: add circumcentres of low-quality faces until none remain.
    # (restarts the scan after every addition, since t.f is mutated in place.)
    idle = false
    while !idle
        idle = true
        for fc in t.f
            c, r = circumcircle_data(fc)
            minlen = minimum(edge_length, neighbours(fc))
            if r / minlen > quality
                add_to_triangulation!(t, fc, c; delaunay=true)
                idle = false
                break
            end
        end
    end

    for i in n+1:length(t.v)
        t.v[i].is_fake = true
    end

    t
end

################################################################
# equidistributed points, wiggling (triangulations.gi:626-720)

"""N P1 points, reasonably well spaced over the sphere, over the field F."""
function equidistributed_p1_points(N::Integer, F::Field)
    Fr = real_field(zero(F)) # P1Point_R3 needs real (not complex-capable) coordinates
    p = P1Point{elem_type(F)}[]
    while length(p) < min(N,10) # add a little randomness
        x = (rand(-10^5:10^5), rand(-10^5:10^5), rand(-10^5:10^5))
        x == (0,0,0) && continue # that would be VERY unlucky
        push!(p, P1Point_R3(Fr.(x)))
    end
    length(p) == N && return p

    t = delaunay_triangulation(p)
    r = pi/2
    while length(t.v) < N
        i = findfirst(fc -> radius(fc) >= r, t.f)
        if i === nothing
            r *= 3/4
            continue
        end
        add_to_triangulation!(t, t.f[i], centre(t.f[i]); delaunay=true)
    end
    pos.(t.v)
end

"""a new triangulation with only the P1 coordinates changed: `movement` is
either a Vector of new vertex positions (fake vertices keep their old
position), or a MoebiusP1Map applied to every vertex, or `nothing` (a plain
topological copy). The Delaunay condition is restored by flipping afterwards."""
function wiggled_triangulation(t::SphereTriangulation, movement)
    nv, ne, nf = length(t.v), length(t.e), length(t.f)
    rv = [TVertex(i) for i in 1:nv]
    re = [TEdge(i) for i in 1:ne]
    rf = [TFace(i) for i in 1:nf]

    for i in 1:nv
        rv[i].neighbour = re[neighbour(t.v[i]).index]
        rv[i].is_fake = is_fake(t.v[i])
    end
    for i in 1:ne
        oe = t.e[i]
        re[i].from = rv[edge_from(oe).index]
        re[i].left = rf[edge_left(oe).index]
        re[i].nextedge = re[nextedge(oe).index]
        re[i].prevopp = re[prevopp(oe).index]
        re[i].opposite = re[opposite(oe).index]
        has_group_element(oe) && (re[i].groupelement = group_element(oe))
    end
    for i in 1:nf
        rf[i].neighbour = re[neighbour(t.f[i]).index]
    end

    r = SphereTriangulation(rv, re, rf)

    if movement isa AbstractVector
        for i in 1:nv
            rv[i].pos = is_fake(rv[i]) ? pos(t.v[i]) : movement[i]
        end
    elseif movement isa MoebiusP1Map
        for i in 1:nv
            rv[i].pos = image(movement, pos(t.v[i]))
        end
    else
        for i in 1:nv
            rv[i].pos = pos(t.v[i])
        end
    end

    fix_delaunay!(r)
    r
end

shallow_copy(t::SphereTriangulation) = wiggled_triangulation(t, nothing)

################################################################
# closest object dispatches (triangulations.gi:706-720)

closest_face(x::TVertex) = edge_left(neighbour(x))
closest_face(x::TEdge) = edge_left(x)
closest_face(x::TFace) = x

closest_faces(x::TVertex) = [edge_left(e) for e in neighbours(x)]
closest_faces(x::TEdge) = [edge_left(x), edge_right(x)]
closest_faces(x::TFace) = [x]

closest_vertex(x::TVertex) = x
closest_vertex(x::TEdge) = edge_from(x)
closest_vertex(x::TFace) = edge_from(neighbour(x))

closest_vertices(x::TVertex) = [x]
closest_vertices(x::TEdge) = [edge_from(x), edge_to(x)]
closest_vertices(x::TFace) = [edge_from(e) for e in neighbours(x)]
