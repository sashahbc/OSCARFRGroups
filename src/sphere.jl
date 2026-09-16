# fundamental groups of punctured spheres
#
# A sphere group is the fundamental group of a punctured sphere, or more
# generally an orbifold whose underlying space is a sphere: generators
# x_1,...,x_n with relation x_{ordering[1]}*...*x_{ordering[n]} = 1, and
# (optionally) x_i^{exponents[i]} = 1.
#
# The group-theoretic machinery here -- a custom Knuth-Bendix rewriting
# system for sphere-group normal forms (with special-cased exceptional
# orbifold configurations), curve intersection numbers (Cohen-Lustig),
# and mapping-class-group / sphere-automorphism-group computation
# (braid-twist factorization) -- is deep, specialized combinatorial group
# theory that GAP's "img" package already implements. Rather than
# reimplementing it, SphereGroup wraps GAP's own SphereGroup type directly,
# the same way FRGroup wraps GAP's FR groups in FRGroups.jl, with a
# Julia-idiomatic API on top.
#
# This requires the "img" GAP package to be loadable (GAP.Packages.load
# ("IMG"), invoked below in __init__ via FRGroups.jl) -- e.g. via a symlink
# at ~/.gap/pkg/img pointing at a checkout of
# https://github.com/gap-packages/img (the compiled kernel extension,
# img.so, is NOT required for sphere groups -- only for the P1-point
# machinery elsewhere in that package).

export SphereGroup, sphere_group, SphereGroupElem, SphereConjClass
export sphere_rank, sphere_ordering, sphere_exponents
export isomorphism_free_group
export as_sphere_group, isomorphism_sphere_group
export peripheral_classes, is_peripheral
export intersection_number, self_intersection_number
export sphere_automorphism_group, epimorphism_to_out
export amalgamated_free_product, embeddings_of_amalgamated_free_product

@attributes mutable struct SphereGroup <: Oscar.GAPGroup
    X::GapObj
    function SphereGroup(G::GapObj)
        @assert GAP.Globals.IsSphereGroup(G)
        return new(G)
    end
end

const SphereGroupElem = Oscar.BasicGAPGroupElem{SphereGroup}
const SphereConjClass = Oscar.GAPGroupConjClass{SphereGroup,SphereGroupElem}

"""
    sphere_group(ordering::Integer, exponents=nothing)
    sphere_group(ordering::AbstractVector{<:Integer}, exponents=nothing)

The fundamental group of a punctured sphere / sphere orbifold: generators
`x_1,...,x_n` with relation `x_{ordering[1]}*...*x_{ordering[n]} = 1`, and
`x_i^{exponents[i]} = 1` for each `i` with `exponents[i] != 0` (`exponents`
defaults to all-0, i.e. all generators of infinite order -- then the group
is isomorphic to a free group of rank `n-1`, see [`isomorphism_free_group`](@ref)).

`ordering` may be a single integer `n`, meaning the ordering `[n,n-1,...,1]`.
"""
function sphere_group(ordering::Union{Integer,AbstractVector{<:Integer}}, exponents::Union{Nothing,AbstractVector{<:Integer}} = nothing)
    o = ordering isa Integer ? GAP.Obj(Int(ordering)) : GapObj(Int.(collect(ordering)); recursive=true)
    G = exponents === nothing ? GAP.Globals.SphereGroup(o) : GAP.Globals.SphereGroup(o, GapObj(Int.(collect(exponents)); recursive=true))
    SphereGroup(G)
end

Base.show(io::IO, g::SphereGroup) = print(io, String(GAP.Globals.ViewString(GapObj(g))))

"""the number of generators of the sphere group g."""
sphere_rank(g::SphereGroup) = Int(GAP.Globals.RankOfSphereGroup(GapObj(g)))

"""the ordering of g, i.e. the permutation of generators whose product is trivial."""
sphere_ordering(g::SphereGroup) = Vector{Int}(GAP.Globals.OrderingOfSphereGroup(GapObj(g)))

"""the exponents of the generators of g (0 meaning infinite order)."""
sphere_exponents(g::SphereGroup) = Vector{Int}(GAP.Globals.ExponentsOfSphereGroup(GapObj(g)))

"""the Euler characteristic of g (a free group, or a sphere group). g is
finite iff this is positive, and virtually abelian iff this is 0."""
Oscar.euler_characteristic(g::Oscar.GAPGroup) = Rational{Int}(GAP.Globals.EulerCharacteristic(GapObj(g)))

"""if g has all exponents infinite, an isomorphism (as a raw GAP group
homomorphism) from g to a free group of rank sphere_rank(g)-1."""
isomorphism_free_group(g::SphereGroup) = GAP.Globals.IsomorphismFreeGroup(GapObj(g))

"""convert the f.p. group G (which must have exactly one relator that's a
permutation of all the generators, the "ordering" relator, and otherwise
only power relators) into a sphere group; returns `nothing` if it can't."""
function as_sphere_group(G::Oscar.GAPGroup)
    newG = GAP.Globals.AsSphereGroup(GapObj(G))
    newG == GAP.Globals.fail && return nothing
    SphereGroup(newG)
end

"""an isomorphism (as a raw GAP group homomorphism) from the f.p. group G to
a sphere group, or `nothing` if G is not (as given) a sphere group."""
function isomorphism_sphere_group(G::Oscar.GAPGroup)
    iso = GAP.Globals.IsomorphismSphereGroup(GapObj(G))
    iso == GAP.Globals.fail && return nothing
    iso
end

################################################################
# conjugacy classes ("multicurves")
################################################################

"""the peripheral conjugacy classes of g, i.e. the conjugacy classes of its generators."""
peripheral_classes(g::SphereGroup) = [conjugacy_class(g,x) for x=gens(g)]

"""whether x represents a peripheral (boundary) loop, i.e. is conjugate to a
power of a single generator."""
is_peripheral(x::SphereGroupElem) = Bool(GAP.Globals.IsPeripheral(GapObj(x)))
is_peripheral(c::SphereConjClass) = Bool(GAP.Globals.IsPeripheral(GapObj(c)))

"""the geometric intersection number of the (multi)curves represented by the
conjugacy classes c,d: the minimal number of intersections they may have,
following Cohen-Lustig."""
intersection_number(c::SphereConjClass, d::SphereConjClass) = Int(GAP.Globals.IntersectionNumber(GapObj(c), GapObj(d)))

"""the self-intersection number of the curve represented by c, i.e. its
intersection number with a small translate of itself."""
self_intersection_number(c::SphereConjClass) = Int(GAP.Globals.SelfIntersectionNumber(GapObj(c)))

################################################################
# automorphism groups, amalgamated free products
################################################################

"""the group (as a raw GAP group) of automorphisms of the sphere group g
that preserve all its peripheral conjugacy classes."""
sphere_automorphism_group(g::SphereGroup) = GAP.Globals.AutomorphismGroup(GapObj(g))

"""an epimorphism (as a raw GAP group homomorphism) from the sphere
automorphism group a (as returned by sphere_automorphism_group) to its
group of outer automorphisms."""
epimorphism_to_out(a::GapObj) = GAP.Globals.EpimorphismToOut(a)

"""the amalgamated free product of the sphere groups g and h, along the
cyclic subgroups generated by x∈g and y∈h."""
function amalgamated_free_product(g::SphereGroup, h::SphereGroup, x::SphereGroupElem, y::SphereGroupElem)
    SphereGroup(GAP.Globals.AmalgamatedFreeProduct(GapObj(g), GapObj(h), GapObj(x), GapObj(y)))
end

"""the two embeddings (as raw GAP group homomorphisms) of the original
factors g,h into amal = amalgamated_free_product(g,h,x,y). Must be called
on amal itself, not on g or h."""
embeddings_of_amalgamated_free_product(amal::SphereGroup) = Vector{GapObj}(GAP.Globals.EmbeddingsOfAmalgamatedFreeProduct(GapObj(amal)))
