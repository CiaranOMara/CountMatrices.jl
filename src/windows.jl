abstract type AbstractGenomicWindow{T} <: GenomicFeatures.AbstractGenomicInterval{T} end

"Window with fixed radius or width."
struct GenomicWindowCentred{R, M, C<:GenomicFeatures.AbstractGenomicInterval} <: AbstractGenomicWindow{M}
	centre::C
	metadata::M
end

function GenomicWindowCentred(centre::C, metadata::M, radius::Integer) where {M, C<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindowCentred{radius, M, C}(centre, metadata)
end

function GenomicWindowCentred(centre::C, radius::Integer) where {C<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindowCentred(centre, nothing, radius)
end

function GenomicWindowCentred{R, M}(centre::C) where {R, M, C<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindowCentred(centre, M(), R)
end

function GenomicWindowCentred{R}(centre::C, mdata) where {R, C<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindowCentred(centre, mdata, R)
end

function GenomicWindowCentred{R}(centre::C) where {R, C<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindowCentred(centre, nothing, R)
end

"Radius of the window."
function radius(::Type{GenomicWindowCentred{R, M, C}}) where {R, M, C}
	return R
end

function radius(window::GenomicWindowCentred{R, M, C}) where {R, M, C}
	return R
end

function radius(::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {R, M, C}
	return R
end

function GenomicFeatures.leftposition(window::GenomicWindowCentred)
	return leftposition(window.centre) - radius(window)
end

function GenomicFeatures.rightposition(window::GenomicWindowCentred)
	return rightposition(window.centre) + radius(window)
end

"Window with dynamic/variable left and right positions."
struct GenomicWindow{M, C<:GenomicFeatures.AbstractGenomicInterval} <: AbstractGenomicWindow{M}
	centre::C
	first::Int64
	last::Int64
	metadata::M
end

function GenomicWindow(centre::C, left::Integer, right::Integer, metadata::M) where {M, C<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindow{M, C}(centre, left, right, metadata)
end

function GenomicWindow(centre::C, left::Integer, right::Integer) where {C<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindow(centre, left, right, nothing)
end

function GenomicWindow(centre::C, interval::I) where {C<:GenomicFeatures.AbstractGenomicInterval, I<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindow(centre, leftposition(interval), rightposition(interval), nothing)
end

function GenomicWindow(centre::C, interval::I, metadata) where {C<:GenomicFeatures.AbstractGenomicInterval, I<:GenomicFeatures.AbstractGenomicInterval}
	return GenomicWindow(centre, leftposition(interval), rightposition(interval), metadata)
end

function GenomicFeatures.groupname(window::AbstractGenomicWindow)
	return groupname(window.centre)
end

function centre(window::AbstractGenomicWindow)
	return window.centre
end

function Base.:(==)(a::AbstractGenomicWindow, b::AbstractGenomicWindow)
	return centre(a) == centre(b) &&
		leftposition(a) == leftposition(b) &&
		rightposition(a) == rightposition(b) &&
		GenomicFeatures.metadata(a) == GenomicFeatures.metadata(b)
end

function Base.empty!(window::AbstractGenomicWindow)
	return empty!(window.metadata)
end

function Base.isempty(n::Number, window::AbstractGenomicWindow)
	return isempty(window.metadata)
end

function Base.push!(window::AbstractGenomicWindow, interval::GenomicFeatures.AbstractGenomicInterval)
	return push!(window.metadata, interval)
end

function Base.append!(window::AbstractGenomicWindow, intervals::AbstractVector{<:GenomicFeatures.AbstractGenomicInterval})
	return append!(window.metadata, intervals)
end
