abstract type Biset end

struct LeftFreeBiset <: Biset
    left::Group
    right::Group
    table::Matrix{Pair{Int,Int}}
end
