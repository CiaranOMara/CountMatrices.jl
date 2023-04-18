function centreposition(interval::GenomicFeatures.AbstractGenomicPosition)
	return leftposition(interval)
end

function centreposition(interval::GenomicFeatures.AbstractGenomicInterval) :: Int
	return leftposition(interval) + round(Int, span(interval) / 2)
end

function centreposition(window::AbstractGenomicWindow)
	return centreposition(centre(window))
end

function relative(n::Number, window::AbstractGenomicWindow)
	return n - centreposition(window)
end

function group!(windows::GenomicFeatures.IterableGenomicCollection{I}, signal;
		skip::Function = (record) -> false,
		preprocess::Function = (record) -> record,
		process::Function = (record) -> record,
	) where I<:AbstractGenomicWindow

	for (window, record) in eachoverlap(windows, signal)
		record = preprocess(record)
		if skip(record)
			continue
		end
		record = process(record)
		push!(window, record)
	end
	return windows
end

function determine_allocation(window::AbstractGenomicWindow{I}) where {T <: Int, I <: GenomicFeatures.AbstractGenomicInterval{T}} # Note: allocation expects coverage intervals.
	return sum(span.(GenomicFeatures.metadata(window)))
end

function determine_allocation(windows::GenomicFeatures.IterableGenomicCollection{<:AbstractGenomicWindow})

	spans = Iterators.Generator(windows) do window
		return span.(GenomicFeatures.metadata(window))
	end |> Iterators.flatten

	return sum(spans)
end

function determine_N(windows::GenomicFeatures.IterableGenomicCollection{<:AbstractGenomicWindow})

	lengths = Iterators.Generator(windows) do window
		return length(GenomicFeatures.metadata(window))
	end |> Iterators.flatten

	return sum(lengths)
end

function relative_xy(windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {I<:GenomicFeatures.AbstractGenomicPosition, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}

	N = determine_N(windows)

	xs = Vector{Int}(undef, N)
	ys = Vector{Float64}(undef, N)

	i = 1

	for window in windows

		for interval in GenomicFeatures.metadata(window)

			xs[i] = relative(leftposition(interval), window)
			ys[i] = GenomicFeatures.metadata(interval)

			i = i + 1
		end

	end

	return xs, ys

end

function relative_xy(windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {I<:GenomicFeatures.AbstractGenomicInterval, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}

	N = determine_allocation(windows) # Note: assuming interval metadata contains coverage.

	xs = Vector{Int}(undef, N)
	ys = Vector{Float64}(undef, N)

	n = 0

	for window in windows # Note: assuming windows are ordered.
		for interval in GenomicFeatures.metadata(window)
			start = relative(leftposition(interval), window) # Note: filtering at the end.
			stop = relative(rightposition(interval), window)

			r = range(start; stop)

			l = length(r)

			lo = n + 1
			hi = n + l

			xs[lo:hi] = r
			ys[lo:hi] .= GenomicFeatures.metadata(interval)

			n = hi
		end

	end

	select = (xs .>= -R) .& (xs .<= R)

	return xs[select], ys[select]

end

function relative_xmin_xmax_y(windows::GenomicFeatures.IterableGenomicCollection{<:GenomicWindowCentred{R, M, C}}) where {I<:GenomicFeatures.AbstractGenomicInterval, R, M<:GenomicFeatures.IterableGenomicCollection{I}, C}

	N = determine_N(windows)

	xmins = Vector{Int}(undef, N)
	xmaxs = Vector{Int}(undef, N)
	ys = Vector{Float64}(undef, N)

	i = 1

	for window in windows

		for interval in GenomicFeatures.metadata(window)

			xmins[i] = relative(leftposition(interval), window)
			xmaxs[i] = relative(rightposition(interval), window)

			ys[i] = GenomicFeatures.metadata(interval)

			i = i + 1
		end

	end

	return xmins, xmaxs, ys

end

function GenomicFeatures.coverage(window::GenomicWindowCentred)
	return GenomicWindowCentred(centre(window), coverage(GenomicFeatures.metadata(window)), radius(window))
end

function zero_nonfinite!(cm::AbstractCountMatrix{T}) where T
	select = isfinite.(cm.m) .== false
	cm.m[select] .= zero(T)
	return cm
end
