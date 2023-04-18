abstract type AbstractCountMatrix{T} end

struct CountMatrix{T} <: AbstractCountMatrix{T}
	m::Matrix{T}
	oa::OffsetArray
end

function CountMatrix(windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {T<:Number, I<:GenomicFeatures.AbstractGenomicInterval{T}, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}

	domain = -R:R
	seqplot_range = 1:length(windows) #Note: The range symbol is occupied by Base.range.

	m = fill(NaN, length(domain), length(seqplot_range)) #Note: row, column. Windows run along columns.
	oa = OffsetArray(m, domain, seqplot_range)

	for (i, window) in enumerate(windows) # Note: assuming windows are ordered.
		for interval in GenomicFeatures.metadata(window)
			start = max(-R, relative(leftposition(interval), window))
			stop = min(R, relative(rightposition(interval), window))
			oa[range(start; stop), i] .= GenomicFeatures.metadata(interval) #Note: assuming interval metadata contains coverage.
		end
	end

	return CountMatrix(m, oa)
end

function CountMatrix(windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {I<:GenomicFeatures.AbstractGenomicInterval{Nothing}, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}
	# new_windows = Iterators.Generator(windows) do window
	# 	return GenomicWindowCentred{R}(centre(window), coverage(GenomicFeatures.metadata(window)))
	# end
	#
	# return CountMatrix(collect(new_windows)) #Note: Generator does not expose type. We also need to collect to avoid Base.Generator{Base.Generatror... ect} loop.

	# iter = Base.Generator(coverage, windows)
	# return CountMatrix(iter) # no method matching CountMatrix(::Base.Generator{Base.Generator{Vector{GenomicWindowCentred{2, Vector{GenomicInterval{Nothing}}, GenomicPosition{Nothing}}}, typeof(coverage)}, typeof(coverage)})

	return CountMatrix(coverage.(windows))

end

struct MemoryMappedCountMatrix{T, I<:IO} <: AbstractCountMatrix{T}
	io::I
	m::Matrix{T}
	oa::OffsetArray
end

function MemoryMappedCountMatrix(io::IO, windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {T<:Number, I<:GenomicFeatures.AbstractGenomicInterval{T}, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}

	#TODO: ensure unassigned are NaNs.

	domain = -R:R
	seqplot_range = 1:length(windows) #Note: The range symbol is occupied by Base.range.

	write(io, length(domain))
	write(io, length(seqplot_range))

	m = Mmap.mmap(io, Matrix{Float64}, (length(domain), length(seqplot_range)))
	oa = OffsetArray(m, domain, seqplot_range)

	for (i, window) in enumerate(windows) # Note: assuming windows are ordered.
		for interval in GenomicFeatures.metadata(window)
			start = max(-R, relative(leftposition(interval), window))
			stop = min(R, relative(rightposition(interval), window))
			oa[range(start; stop), i] .= GenomicFeatures.metadata(interval) #Note: assuming interval metadata contains coverage.
		end
	end

	Mmap.sync!(m)

	return MemoryMappedCountMatrix(io, m, oa)
end

function MemoryMappedCountMatrix(io::IO, windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {I<:GenomicFeatures.AbstractGenomicInterval{Nothing}, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}
	return MemoryMappedCountMatrix(io, coverage.(windows))
end

function MemoryMappedCountMatrix(windows::GenomicFeatures.IterableGenomicCollection)
	path = tempname(@get_scratch!("countmatrix"))
	return MemoryMappedCountMatrix(open(path, "w+"), windows) #TODO: is a record of path needed for saving and copying work?
end

struct SparseCountMatrix{Tv} <: AbstractCountMatrix{Tv}
	m::SparseMatrixCSC{Tv,<:Integer}
	oa::OffsetArray
end

function SparseCountMatrix(columns, rows, values::AbstractVector{T}) where T
	domain = minimum(columns):maximum(columns)
	nrows = length(unique(rows))

	m = sparse(columns, rows, values)
	oa = OffsetArray(m, domain, 1:nrows)

	return SparseCountMatrix{T}(m, oa)
end

function SparseCountMatrix(windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {T<:Number, I<:GenomicFeatures.AbstractGenomicInterval{T}, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}

	#TODO: ensure unassigned are NaNs.

	N = determine_allocation(windows) # Note: assuming interval metadata contains coverage.

	xs = Vector{Int}(undef, N)
	ys = Vector{Int}(undef, N)
	zs = Vector{T}(undef, N)

	n = 0
	for (j, window) in enumerate(windows) # Note: assuming windows are ordered.
		for interval in GenomicFeatures.metadata(window)
			start = relative(leftposition(interval), window) # Note: filtering at the end.
			stop = relative(rightposition(interval), window)

			r = range(start; stop)

			l = length(r)

			lo = n + 1
			hi = n + l

			xs[lo:hi] = r
			ys[lo:hi] .= j
			zs[lo:hi] .= GenomicFeatures.metadata(interval)

			n = hi
		end
	end

	# select = Base.between.(xs, -R, R)
	select = (xs .>= -R) .& (xs .<= R)

	domain = -R:R
	seqplot_range = 1:length(windows) #Note: The range symbol is occupied by Base.range.

	m = sparse(xs[select] .+ (R +1) , ys[select], zs[select]) #TODO: determine index selection returns a copy.
	oa = OffsetArray(m, domain, seqplot_range)

	return SparseCountMatrix{T}(m, oa)
end

function SparseCountMatrix(windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {I<:GenomicFeatures.AbstractGenomicInterval{Nothing}, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}
	return SparseCountMatrix(coverage.(windows))
end

function domain(cm::AbstractCountMatrix)
	return axes(cm.oa, 1)
end

function Base.range(cm::AbstractCountMatrix)
	return axes(cm.oa, 2)
end

function Base.eachcol(cm::AbstractCountMatrix)
	return eachcol(cm.m) #Note: bypassing OffsetArray.
end

function Base.eachrow(cm::AbstractCountMatrix)
	return eachrow(cm.oa)
end

#TODO: dynamically construct array method passthroughs with meta programming.
function Base.extrema(cm::AbstractCountMatrix)
	return extrema(cm.m)
end
