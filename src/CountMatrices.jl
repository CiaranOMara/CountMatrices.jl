module CountMatrices

using Mmap
using Scratch
using SparseArrays

using GenomicFeatures
using OffsetArrays


export CountMatrix, SparseCountMatrix, MemoryMappedCountMatrix, GenomicWindow, GenomicWindowCentred, centre

Broadcast.broadcastable(x::GenomicFeatures.IntervalTrees.AbstractInterval) = Ref(x)

include("windows.jl")

include("matrices.jl")

include("summaries.jl")

include("utils.jl")

end
