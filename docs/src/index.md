```@meta
CurrentModule = CountMatrices
```

# CountMatrices

Documentation for [CountMatrices](https://github.com/CiaranOMara/CountMatrices.jl).

```@example quickstart
using CountMatrices

using FileIO
using GenomicFeatures # v3
using ColorSchemes

# Define GenomicIntervals.
intervals = [
	GenomicPosition("test1", 1, 1),
	GenomicPosition("test2", 2, 2),
	GenomicPosition("test3", 1, 3),
	GenomicPosition("test3", 3, 3),
]

# Define various genomic windows.
windows = [
	GenomicWindowCentred(GenomicPosition("test1", 2), Vector{eltype(intervals)}(), 1),
	GenomicWindowCentred{1}(GenomicPosition("test2", 2), Vector{eltype(intervals)}()),
	GenomicWindowCentred{1, Vector{eltype(intervals)}}(GenomicPosition("test3", 2)),
]

# Collect intervals overlapping windows.
for (window, interval) in eachoverlap(windows, intervals)
	push!(window, interval)
end

# Sort windows to place the most significant at the top of the data matrix.
sort!(windows; by=window->sum(GenomicFeatures.volume, GenomicFeatures.metadata(window)), rev = true)

# Splice windows into count matrix.
M = CountMatrix(windows)

# Sort windows to place the most significant at the top of the data matrix.
# M.m .= sortslices(M.m, by = sum, dims = 2, rev = true)
```

```@example quickstart
palette = ColorSchemes.balance
```

```@example quickstart
# Get RGBA matrix.
path, io, c = CountMatrices.draw(M; palette, rangescale = :extrema, scale = log10)

# Show RGBA matrix.
c
```

# API

```@index
```

```@autodocs
Modules = [CountMatrices]
```
