module FRGroups

using GAP, Oscar, Buchi, NonlinearSolve

export FRGroup, fr_group, mealy_machine

function __init__()
    GAP.Packages.load("fr")

    GAP.Globals.Read(GapObj("$(@__DIR__)/../gap/buchi.gd")) # add new julia automata to FR
    GAP.Globals.Read(GapObj("$(@__DIR__)/../gap/buchi.gi"))

    # needed for SphereGroup (sphere.jl); see that file for how to make the
    # "img" GAP package discoverable if this fails to load.
    if !GAP.Packages.load("IMG")
        @warn "The GAP package \"IMG\" could not be loaded; SphereGroup and related functions will not work. See src/sphere.jl for setup instructions."
    end
end

# ... add a bunch of Julia methods that allow direct creation, in Julia syntax, of FRGroups etc.

# wrapper around FR GAP groups
@attributes mutable struct FRGroup <: Oscar.GAPGroup
   X::GapObj

   function FRGroup(G::GapObj)
     @assert GAP.Globals.IsFRGroup(G)
     return new(G)
   end
end

fr_group(G::GapObj) = FRGroup(G)

const FRGroupElem = Oscar.BasicGAPGroupElem{FRGroup}

include("mealy.jl")

include("bisets.jl")

include("sphere.jl")

include("p1points.jl")

include("triangulations.jl")

include("spider.jl")

include("markedsphere.jl")

include("hurwitz.jl")

include("thurston.jl")

include("examples.jl")

end # module FRGroups
