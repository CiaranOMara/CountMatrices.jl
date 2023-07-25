# CountMatrices.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://CiaranOMara.github.io/CountMatrices.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CiaranOMara.github.io/CountMatrices.jl/dev/)
[![Build Status](https://github.com/CiaranOMara/CountMatrices.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/CiaranOMara/CountMatrices.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/CiaranOMara/CountMatrices.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/CiaranOMara/CountMatrices.jl)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/C/CountMatrices.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/report.html)

The `CountMatrices` package handles the gathering and splicing of metadata into a count matrix.

Example workflow:

```julia
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

# Get RGBA matrix.
path, io, c = draw(M; palette = ColorSchemes.balance, rangescale = :extrema, scale = log10)

# Save colour matrix to file using FileIO.
save("heatmap.png", c)
```
