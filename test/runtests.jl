using CountMatrices
using Test

using GenomicFeatures

@testset "CountMatrices.jl" begin

    @testset "Centre" begin
		@test CountMatrices.centreposition(GenomicPosition("chr1", 1)) == 1
		@test CountMatrices.centreposition(GenomicInterval("chr1", 1, 1)) == 1
		@test CountMatrices.centreposition(GenomicInterval("chr1", 1, 9)) == 5
		@test CountMatrices.centreposition(GenomicInterval("chr1", 1, 10)) == 6
	end # Centre

	@testset "GenomicWindow" begin

		gcentre = GenomicPosition("chr1", 10)
		window1 = GenomicWindowCentred(gcentre, Vector{GenomicInterval{Int}}(), 4)
		window2 = GenomicWindowCentred(gcentre, Vector{GenomicInterval{Int}}(), 4)

		intervals = GenomicInterval.("chr1", [8, 11, 19], [9, 12, 21], '.', [7, 8, 9])

		# Test push.
		foreach(intervals) do interval
			push!(window1.metadata, interval)
			push!(window2, interval)
		end

		@test length(window1.metadata) == length(intervals)
		@test length(window2.metadata) == length(intervals)
		@test window1 == window2

		empty!(window1)
		empty!(window2)

		@test length(window1.metadata) == 0
		@test length(window2.metadata) == 0

		# Check that intervals still exist. #TODO: append to view or vector of references.
		@test length(intervals) > 0

		# Test append!,
		append!(window1.metadata, intervals)
		append!(window2, intervals)

		@test window1 == window2

		#Check shorthand constructor.
		@test GenomicWindowCentred{4, Vector{GenomicInterval{Int}}}(gcentre) == GenomicWindowCentred(gcentre, Vector{GenomicInterval{Int}}(), 4)

		gcentres =  GenomicPosition.("chr1", [10, 20, 30])
		@test GenomicWindowCentred{4, Vector{GenomicInterval{Nothing}}}.(gcentres) == GenomicWindowCentred.(gcentres, [Vector{GenomicInterval{Int}}()], 4)

	end

	@testset "IterableGenomicWindow" begin

		window = GenomicWindowCentred(GenomicPosition("chr1", 10), Vector{GenomicInterval{Int}}(), 4)

		# Note: generators do not provide eltype.
		# iter1 = Base.Generator(x->x, GenomicIntervalCollection([window]))
		# iter2 = Base.Generator(x->x, [window])
		# @test eltype(iter1) == eltype(iter2) == typeof(window)

		@test typeof([window]) <: GenomicFeatures.IterableGenomicCollection
		@test typeof(Base.Generator(x->x, [window])) <: GenomicFeatures.IterableGenomicCollection
		@test typeof(Base.Generator(typeof(window), [window])) <: GenomicFeatures.IterableGenomicCollection
		@test typeof(GenomicIntervalCollection([window])) <: GenomicFeatures.IterableGenomicCollection
		@test typeof(Base.Generator(x->x, GenomicIntervalCollection([window]))) <: GenomicFeatures.IterableGenomicCollection
		@test typeof(Base.Generator(typeof(window), GenomicIntervalCollection([window]))) <: GenomicFeatures.IterableGenomicCollection
	end

	@testset "CountMatrix" begin

		gcentres =  GenomicPosition.(string.("chr", [1, 2, 3]), 5)

		intervals =[
			[
				GenomicInterval("chr1", 1, 9),
			], [
				GenomicInterval("chr2", 1, 5),
				GenomicInterval("chr2", 6, 9),
			], [
				GenomicInterval("chr3", 1, 3),
				GenomicInterval("chr3", 4, 6),
				GenomicInterval("chr3", 7, 9),
			]
		]

		expected_rows = [[1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1]]
		expected_cols = [ones(length(-2:2)), ones(length(-2:2)), ones(length(-2:2))]

		windows = GenomicWindowCentred{2}.(gcentres, intervals)

		scm = SparseCountMatrix(windows)
		cm = CountMatrix(windows)
		mcm = MemoryMappedCountMatrix(windows)

		@test CountMatrices.domain(cm) == -2:2
		@test CountMatrices.domain(scm) == -2:2
		@test CountMatrices.domain(mcm) == -2:2

		@test eachcol(scm) |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachcol(cm)  |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachcol(mcm) |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachrow(scm) |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect
		@test eachrow(cm)  |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect
		@test eachrow(mcm) |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect


		intervals =[
			[
				GenomicInterval("chr1", 1, 9),
			], [
				GenomicInterval("chr2", 2, 5),
				GenomicInterval("chr2", 6, 8),
			], [
				GenomicInterval("chr3", 3, 3),
				GenomicInterval("chr3", 4, 6),
				GenomicInterval("chr3", 7, 7),
			]
		]

		expected_rows = [[1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1]]
		expected_cols = [ones(length(-2:2)), ones(length(-2:2)), ones(length(-2:2))]

		windows = GenomicWindowCentred{2}.(gcentres, intervals)

		scm = SparseCountMatrix(windows)
		cm = CountMatrix(windows)
		mcm = MemoryMappedCountMatrix(windows)

		@test eachcol(scm) |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachcol(cm)  |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachcol(mcm) |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachrow(scm) |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect
		@test eachrow(cm)  |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect
		@test eachrow(mcm) |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect

		intervals =[
			[
				GenomicInterval("chr1", 1, 9, '.', 1),
			], [
				GenomicInterval("chr2", 1, 9, '.', 2),
			], [
				GenomicInterval("chr3", 1, 9, '.', 3),
			]
		]

		expected_rows = [[1,2,3], [1,2,3], [1,2,3], [1,2,3], [1,2,3]]
		expected_cols = [fill(1, 5), fill(2, 5), fill(3, 5)]

		windows = GenomicWindowCentred{2}.(gcentres, intervals)

		scm = SparseCountMatrix(windows)
		cm = CountMatrix(windows)
		mcm = MemoryMappedCountMatrix(windows)

		@test eachcol(scm) |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachcol(cm)  |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachcol(mcm) |> Iterators.Flatten |> collect == expected_cols |> Iterators.Flatten |> collect
		@test eachrow(scm) |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect
		@test eachrow(cm)  |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect
		@test eachrow(mcm) |> Iterators.Flatten |> collect == expected_rows |> Iterators.Flatten |> collect

	end

end
