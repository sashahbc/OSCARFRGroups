# fundamental groups of punctured spheres

abstract type SphereGroup <: Group end

struct RWSSphereGroup{T} <: SphereGroup
    rank::Int
    orders::Vector{Int}
    rws::Any
end
