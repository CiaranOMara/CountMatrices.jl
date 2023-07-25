using ColorSchemes
using Colors
# using ImageTransformations #Note: supplies zero(::ColorTypes.RGB{Float64}). #TODO: check whether ColorVectors.jl supplies required functionality.

using ImageTransformations

function noscaling() end

function getrangescale(x, rangemode::Symbol)
	# based on https://github.com/JuliaGraphics/ColorSchemes.jl/blob/24a51534f875428c5d722356edfb8d8b9bb3800b/src/ColorSchemes.jl#L274-L285

	itr = Iterators.filter(isfinite, x)

	if rangemode == :clamp
		return ColorSchemes.defaultrange(itr) .|> float
	end

	if rangemode == :extrema
		return extrema(itr) .|> float
	end

	if rangemode == :centered
		return (-1, 1) .* maximum(abs, itr) .|> float
	end

	throw(ArgumentError("rangescale :$rangemode not supported, should be :clamp, :extrema, :centered or tuple (minVal, maxVal)"))

end

function getrangescale(x, rangemode::NTuple{2,Float64}) # Note: Float64 is required by ColorSchemes.
	return rangemode
end

function colourise(m::AbstractMatrix{T}; palette::ColorScheme, rangescale, missing_color::Colorant) where T <: Number

	# Determine colour element type.
	el = promote_type(eltype(palette.colors), typeof(missing_color))

	# Setup scratch file.
	(path_colorised, io_colorised) = mktemp(@get_scratch!("countmatrix-colorised"))

	# Initialise Colour matrix.
	# c = Array{el}(undef, size(m))
	c = mmap(io_colorised, Matrix{el}, size(m))

	# Add colour to matrix column by column.
	for (col_c, col_m) in zip(eachcol(c), eachcol(m)) # Note: col is a view.

		#
		select = isfinite.(col_m)

		# Colorise column.
		col_c[select] .= get(palette, col_m[select], rangescale) # note: converts RGB to RGBA on assignment.

		# Add colour for missing values.
		col_c[select .== false] .= missing_color #TODO: Use min colour for -Inf and NaN, and max colour for +Inf.
	end

	return path_colorised, io_colorised, c
end

function resize(c::AbstractMatrix{T}, dims::Dims{2}) where T <: Color
	# Setup scratch file.
	(path_resized, io_resized) = mktemp(@get_scratch!("countmatrix-colorised-resized"))

	# Initialise matrix with final dimensions.
	resized = mmap(io_resized, Matrix{T}, reverse(dims))

	# Resize image.
	imresize!(resized, c) # TODO: facilitate selection of resizing algorithm.

	return path_resized, io_resized, resized
end

function draw(M::AbstractMatrix; palette::ColorScheme = ColorSchemes.batlow, rangescale = :extrema, scale::Function = noscaling, missing_color::Colorant = first(palette), dims::Tuple = tuple(), kwargs...)

	# Scale matrix values.
	if nameof(scale) == :noscaling
		m = M
	else
		m = scale.(M)
	end

	# Apply range scale to matrix not isolated columns.
	rangescale = getrangescale(m, rangescale)

	# Colourise matrix.
	path_colorised, io_colorised, c = colourise(m; palette, rangescale, missing_color)

	if isempty(dims)
		return path_colorised, io_colorised, c
	end

	# Resize matrix.
	path_resized, io_resized, resized = resize(c, dims)

	close(io_colorised)
	rm(path_colorised)

	return path_resized, io_resized, resized
end

function draw(SM::SparseArrays.AbstractSparseMatrixCSC; kwargs...)
	return draw(Matrix(SM); kwargs...)
end

function draw(CM::CountMatrices.AbstractCountMatrix; kwargs...)
	return draw(transpose(CM.m); kwargs...) # Swapping rows and columns for image.
end
