@testset "Constructors" begin

    """
        make_z2_mt()

    Multiplication table for Z2 with basis:
    1 = vacuum, 2 = nontrivial element, and 2⊗2 = 1.
    """
    function make_z2_mt()
        mt = zeros(Int, 2, 2, 2)
        mt[1,1,1] = 1
        mt[1,2,2] = 1
        mt[2,1,2] = 1
        mt[2,2,1] = 1
        mt
    end

    """
        make_z3_mt()

    Multiplication table for Z3 with basis:
    1↔0, 2↔1, 3↔2 mod 3.
    """
    function make_z3_mt()
        mt = zeros(Int, 3, 3, 3)
        for i in 0:2, j in 0:2
            k = mod(i + j, 3)
            mt[i+1, j+1, k+1] = 1
        end
        mt
    end

    #Core validation:

        @testset "fusion_ring: valid construction" begin
        mt = make_z2_mt()
        r = fusion_ring(mt; labels = ["0", "1"], names = ["Z2"])

        check_equal(rank(r), 2,
            "fusion_ring did not construct a rank-2 ring from a valid Z2 multiplication table")
        check_equal(size(multiplication_table(r)), (2,2,2),
            "fusion_ring did not preserve the multiplication-table dimensions for a valid Z2 table")
        check_equal(labels(r), ["0", "1"],
            "fusion_ring did not preserve labels for a valid Z2 table")
        check_equal(names(r), ["Z2"],
            "fusion_ring did not preserve names for a valid Z2 table")
        check_true(is_commutative(r),
            "fusion_ring constructed a valid Z2 table but the result was not detected as commutative")
        check_true(is_group_ring(r),
            "fusion_ring constructed a valid Z2 table but the result was not detected as a group ring")
    end

    @testset "fusion_ring: default labels are created" begin
        mt = make_z2_mt()
        r = fusion_ring(mt)

        check_equal(rank(r), 2,
            "fusion_ring without explicit labels did not construct a rank-2 ring")
        check_equal(length(labels(r)), 2,
            "fusion_ring without explicit labels did not create exactly 2 labels")
    end

    @testset "fusion_ring: skip_check=true allows raw construction" begin
        mt_bad = zeros(Int, 2, 2, 2)
        # nonsense table, but construction should succeed with skip_check=true
        r = fusion_ring(mt_bad; skip_check = true)

        check_equal(rank(r), 2,
            "fusion_ring(...; skip_check=true) did not construct a ring object from a raw 2×2×2 tensor")
        check_equal(size(multiplication_table(r)), (2,2,2),
            "fusion_ring(...; skip_check=true) did not preserve multiplication-table dimensions")
    end

    @testset "fusion_ring: rejects negative structure constants" begin
        mt = make_z2_mt()
        mt[2,2,1] = -1

        check_throws(
            () -> fusion_ring(mt; labels = ["0", "1"]),
            "fusion_ring accepted a multiplication table with a negative structure constant"
        )
    end

     @testset "fusion_ring: rejects non-integer structure constants" begin
        mt = Array{Float64}(undef, 2, 2, 2)
        fill!(mt, 0.0)
        mt[1,1,1] = 1.0
        mt[1,2,2] = 1.0
        mt[2,1,2] = 1.0
        mt[2,2,1] = 1.0

        check_throws(
            () -> fusion_ring(mt; labels = ["0", "1"]),
            "fusion_ring accepted a multiplication table with non-integer structure constants"
        )
    end

    @testset "fusion_ring: rejects non-cubic tensors" begin
        mt = zeros(Int, 2, 2, 3)

        check_throws(
            () -> fusion_ring(mt),
            "fusion_ring accepted a multiplication table whose tensor dimensions were not all equal"
        )
    end

    @testset "fusion_ring: rejects bad unit" begin
        mt = make_z2_mt()
        mt[1,2,2] = 0

        check_throws(
            () -> fusion_ring(mt; labels = ["0", "1"]),
            "fusion_ring accepted a multiplication table whose first basis element was not a unit"
        )
    end

      @testset "fusion_ring: rejects bad inverse condition" begin
        mt = make_z2_mt()
        # make both 1 and 2 appear as "duals" of 2 by forcing extra contribution to c=1
        mt[2,2,2] = 1

        check_throws(
            () -> fusion_ring(mt; labels = ["0", "1"]),
            "fusion_ring accepted a multiplication table violating the unique inverse condition"
        )
    end

    @testset "fusion_ring: rejects non-associative tables" begin
        mt = zeros(Int, 2, 2, 2)
        mt[1,1,1] = 1
        mt[1,2,2] = 1
        mt[2,1,2] = 1
        mt[2,2,2] = 1   # x⊗x = x, this breaks inverse condition / associativity expectations

        check_throws(
            () -> fusion_ring(mt; labels = ["0", "1"]),
            "fusion_ring accepted a non-fusion multiplication table that should fail validation"
        )
    end

    @testset "fusion_ring: rejects incorrect label length" begin
        mt = make_z2_mt()

        check_throws(
            () -> fusion_ring(mt; labels = ["0"]),
            "fusion_ring accepted labels whose length did not equal the rank"
        )
    end

