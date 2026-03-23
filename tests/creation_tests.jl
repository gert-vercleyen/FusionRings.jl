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

    #Zn fusion rings:

      @testset "zn_fusion_ring" begin
        z1 = zn_fusion_ring(1)
        z2 = zn_fusion_ring(2)
        z3 = zn_fusion_ring(3)
        z4 = zn_fusion_ring(4)

        @testset "basic ranks" begin
            check_equal(rank(z1), 1, "rank(zn_fusion_ring(1)) was not 1")
            check_equal(rank(z2), 2, "rank(zn_fusion_ring(2)) was not 2")
            check_equal(rank(z3), 3, "rank(zn_fusion_ring(3)) was not 3")
            check_equal(rank(z4), 4, "rank(zn_fusion_ring(4)) was not 4")
        end

        @testset "labels" begin
            check_equal(labels(z1), ["0"], "labels(zn_fusion_ring(1)) were incorrect")
            check_equal(labels(z2), ["0", "1"], "labels(zn_fusion_ring(2)) were incorrect")
            check_equal(labels(z3), ["0", "1", "2"], "labels(zn_fusion_ring(3)) were incorrect")
            check_equal(labels(z4), ["0", "1", "2", "3"], "labels(zn_fusion_ring(4)) were incorrect")
        end

        @testset "table sizes" begin
            check_equal(size(multiplication_table(z1)), (1,1,1),
                "multiplication table for zn_fusion_ring(1) had wrong size")
            check_equal(size(multiplication_table(z2)), (2,2,2),
                "multiplication table for zn_fusion_ring(2) had wrong size")
            check_equal(size(multiplication_table(z3)), (3,3,3),
                "multiplication table for zn_fusion_ring(3) had wrong size")
            check_equal(size(multiplication_table(z4)), (4,4,4),
                "multiplication table for zn_fusion_ring(4) had wrong size")
        end

        @testset "basic properties" begin
            check_true(is_commutative(z1), "zn_fusion_ring(1) was not commutative")
            check_true(is_commutative(z2), "zn_fusion_ring(2) was not commutative")
            check_true(is_commutative(z3), "zn_fusion_ring(3) was not commutative")
            check_true(is_commutative(z4), "zn_fusion_ring(4) was not commutative")

            check_true(is_group_ring(z1), "zn_fusion_ring(1) was not detected as a group ring")
            check_true(is_group_ring(z2), "zn_fusion_ring(2) was not detected as a group ring")
            check_true(is_group_ring(z3), "zn_fusion_ring(3) was not detected as a group ring")
            check_true(is_group_ring(z4), "zn_fusion_ring(4) was not detected as a group ring")
        end

        @testset "selected products" begin
            # In Z3: indices 1,2,3 correspond to 0,1,2 mod 3
            check_equal(fusion_product(z3, 1, 1), Dict(1 => 1),
                "vacuum × vacuum was not vacuum in zn_fusion_ring(3)")
            check_equal(fusion_product(z3, 2, 2), Dict(3 => 1),
                "index 2 × index 2 in zn_fusion_ring(3) was not index 3")
            check_equal(fusion_product(z3, 2, 3), Dict(1 => 1),
                "index 2 × index 3 in zn_fusion_ring(3) was not vacuum")

            # In Z4: 2+2 = 0 mod 4? careful with indexing:
            # indices 1,2,3,4 ↔ 0,1,2,3
            # so 3×3 ↔ 2+2 = 0, i.e. vacuum
            check_equal(fusion_product(z4, 3, 3), Dict(1 => 1),
                "index 3 × index 3 in zn_fusion_ring(4) was not vacuum")
        end
    end

    #su2k_fusion_ring:

    @testset "su2k_fusion_ring" begin
        r1 = su2k_fusion_ring(1)
        r2 = su2k_fusion_ring(2)
        r3 = su2k_fusion_ring(3)

        @testset "basic ranks and sizes" begin
            check_equal(rank(r1), 2, "rank(su2k_fusion_ring(1)) was not 2")
            check_equal(rank(r2), 3, "rank(su2k_fusion_ring(2)) was not 3")
            check_equal(rank(r3), 4, "rank(su2k_fusion_ring(3)) was not 4")

            check_equal(size(multiplication_table(r1)), (2,2,2),
                "multiplication table for su2k_fusion_ring(1) had wrong size")
            check_equal(size(multiplication_table(r2)), (3,3,3),
                "multiplication table for su2k_fusion_ring(2) had wrong size")
            check_equal(size(multiplication_table(r3)), (4,4,4),
                "multiplication table for su2k_fusion_ring(3) had wrong size")
        end

        @testset "labels" begin
            check_equal(labels(r1), ["0", "1"],
                "labels(su2k_fusion_ring(1)) were incorrect")
            check_equal(labels(r2), ["0", "1", "2"],
                "labels(su2k_fusion_ring(2)) were incorrect")
            check_equal(labels(r3), ["0", "1", "2", "3"],
                "labels(su2k_fusion_ring(3)) were incorrect")
        end

        @testset "basic properties" begin
            check_true(is_commutative(r1), "su2k_fusion_ring(1) was not commutative")
            check_true(is_commutative(r2), "su2k_fusion_ring(2) was not commutative")
            check_true(is_commutative(r3), "su2k_fusion_ring(3) was not commutative")
        end

        @testset "selected low-k products" begin
            # SU(2)_1 behaves like Z2
            check_equal(fusion_product(r1, 2, 2), Dict(1 => 1),
                "nontrivial simple squared in su2k_fusion_ring(1) was not vacuum")

            # SU(2)_2 has labels 0,1,2.
            # Fusion rule: 1⊗1 = 0 + 2  (indices 2⊗2 = 1 + 3)
            check_equal(fusion_product(r2, 2, 2), Dict(1 => 1, 3 => 1),
                "index 2 × index 2 in su2k_fusion_ring(2) was not vacuum + top object")
        end
    end

    




