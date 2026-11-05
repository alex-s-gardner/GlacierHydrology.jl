using Test
using GlobalGlacierAnalysis
import GlobalGlacierAnalysis as GGA
import GlacierVolumeChange as GVC

@testset "Array Operations" begin
    # `validrange` returns one range per dimension, as a Tuple -- so a 1D input yields a 1-tuple
    # whose only element is the range, not a (first, last) pair.
    @testset "validrange 1D" begin
        # Standard case with valid range
        v = [false, true, true, false, true, false]
        range_result = GGA.validrange(v)
        @test length(range_result) == 1
        @test range_result[1] == 2:5  # first true at 2, last at 5

        # All false -> empty range, not an error
        range_empty = GGA.validrange(falses(10))
        @test isempty(range_empty[1])

        # All true
        @test GGA.validrange(trues(10))[1] == 1:10

        # Single true at start
        @test GGA.validrange([true, false, false])[1] == 1:1

        # Single true at end
        @test GGA.validrange([false, false, true])[1] == 3:3
    end

    @testset "validrange 2D" begin
        # Rectangle of true values embedded in false matrix
        mat = falses(10, 10)
        mat[3:7, 4:8] .= true

        range_2d = GGA.validrange(mat)
        @test length(range_2d) == 2
        @test range_2d[1] == 3:7  # rows
        @test range_2d[2] == 4:8  # columns

        # the returned ranges select exactly the true block back out
        @test all(mat[range_2d...])
        @test sum(mat[range_2d...]) == sum(mat)

        # all false in 2D
        empty_2d = GGA.validrange(falses(4, 4))
        @test isempty(empty_2d[1]) && isempty(empty_2d[2])
    end

    # `validgaps` returns a same-length mask, true only at invalid points *between* the first and
    # last valid point. Leading/trailing invalid runs are boundary, not gaps -- so the mask comes
    # back all-false rather than empty.
    @testset "validgaps" begin
        # Gap in middle
        gaps = GVC.validgaps(BitVector([true, false, false, true, true]))
        @test gaps == BitVector([false, true, true, false, false])
        @test count(gaps) == 2

        # No gaps (all true)
        gaps_none = GVC.validgaps(BitVector([true, true, true, true]))
        @test length(gaps_none) == 4
        @test !any(gaps_none)

        # Leading false values are boundary, not a gap
        @test !any(GVC.validgaps(BitVector([false, false, true, true, true])))

        # Trailing false values are boundary, not a gap
        @test !any(GVC.validgaps(BitVector([true, true, true, false, false])))

        # Multiple gaps
        gaps_multi = GVC.validgaps(BitVector([true, false, true, false, false, true]))
        @test gaps_multi == BitVector([false, true, false, true, true, false])
        @test count(gaps_multi) == 3

        # No valid points at all -> nothing is a gap
        @test !any(GVC.validgaps(falses(5)))
    end

end
