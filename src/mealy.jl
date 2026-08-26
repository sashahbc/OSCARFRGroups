function mealy_machine(transitions::Matrix{Int},output::Vector{PermGroupElem})
    @assert length(output)==size(transitions,1)
    
    GAP.Globals.MealyMachine(GapObj([transitions[i,:] for i=1:size(transitions,1)],recursive=true),GapObj(output,recursive=true))
    # encapsulate it in a Julia object
end

# m = FRGroups.mealy_machine([3 2;3 1;3 2],[cperm([1,2]),cperm(),cperm()])

#function output(m::
