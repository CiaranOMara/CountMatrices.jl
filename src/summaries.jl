function summary(f::Function, m::AbstractMatrix; zero_nonfinite::Bool = false, skip_nonfinite::Bool = false, skip_zeros::Bool = false)

	if !zero_nonfinite && skip_nonfinite
		error("Mutually exclusive arguments given.")
	end

	itr = eachrow(m)

	if zero_nonfinite
		itr = Base.Generator(itr) do row
			new_row = copy(row) # Note: copying to preserve matrix.
			select = isfinite.(new_row) == false
			new_row[select] .= 0.0
			return new_row
		end
	end

	if skip_nonfinite
		itr = Base.Generator(itr) do row
			return Iterators.filter(isfinite, row) # Note: lazy iterator does not provide length for broadcast.
		end
	end

	if skip_zeros
		itr = Base.Generator(itr) do row
			return Iterators.filter(!iszero, row) # Note: lazy iterator does not provide length for broadcast.
		end
	end

	# Collect rows to make them broadcastable.
	itr = Base.Generator(itr) do row
		return collect(row)
	end

	return f.(itr)
end

function summary(f::Function, cm::AbstractCountMatrix; kwargs...)
	return summary(f, cm.m; kwargs...)
end

function summary(v, windows::GenomicFeatures.IterableGenomicCollection{<:AbstractGenomicWindow}; matrixtype::AbstractCountMatrix = CountMatrix, kwargs...)
	cm = matrixtype(windows)
	return summary(v, cm; kwargs...)
end
