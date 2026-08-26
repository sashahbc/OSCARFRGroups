module FRGroups

using GAP, Oscar, Buchi

export FRGroup, fr_group, mealy_machine

function __init__()
    GAP.Packages.load("fr")

    GAP.Globals.Read(GapObj("$(@__DIR__)/../gap/buchi.gd")) # add new julia automata to FR
    GAP.Globals.Read(GapObj("$(@__DIR__)/../gap/buchi.gi"))

    global ℂ, ℂ⁽ᶻ⁾, 𝓏
    ℂ = QQBarField()
    ℂ⁽ᶻ⁾, z = rational_function_field(ℂ,:z)
    𝓏 = p1_map(z)
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
