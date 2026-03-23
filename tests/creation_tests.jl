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

    #psu2k_fusion_ring:

    
    @testset "psu2k_fusion_ring" begin
        r2 = psu2k_fusion_ring(2)
        r4 = psu2k_fusion_ring(4)
        r6 = psu2k_fusion_ring(6)

        @testset "basic ranks and sizes" begin
            check_equal(rank(r2), 2, "rank(psu2k_fusion_ring(2)) was not 2")
            check_equal(rank(r4), 3, "rank(psu2k_fusion_ring(4)) was not 3")
            check_equal(rank(r6), 4, "rank(psu2k_fusion_ring(6)) was not 4")

            check_equal(size(multiplication_table(r2)), (2,2,2),
                "multiplication table for psu2k_fusion_ring(2) had wrong size")
            check_equal(size(multiplication_table(r4)), (3,3,3),
                "multiplication table for psu2k_fusion_ring(4) had wrong size")
            check_equal(size(multiplication_table(r6)), (4,4,4),
                "multiplication table for psu2k_fusion_ring(6) had wrong size")
        end

        @testset "basic properties" begin
            check_true(is_commutative(r2), "psu2k_fusion_ring(2) was not commutative")
            check_true(is_commutative(r4), "psu2k_fusion_ring(4) was not commutative")
            check_true(is_commutative(r6), "psu2k_fusion_ring(6) was not commutative")
        end

        @testset "selected product" begin
            # PSU(2)_2 should also be rank 2 and behave like Z2
            check_equal(fusion_product(r2, 2, 2), Dict(1 => 1),
                "nontrivial simple squared in psu2k_fusion_ring(2) was not vacuum")
        end
    end


    #son2_fusion_ring:

    
    @testset "son2_fusion_ring / metaplectic_fusion_ring" begin
        odd = son2_fusion_ring(5)
        even2 = son2_fusion_ring(6)
        even4 = son2_fusion_ring(8)
        meta = metaplectic_fusion_ring(5)

        @testset "basic sizes" begin
            check_true(rank(odd) > 0, "rank(son2_fusion_ring(5)) was not positive")
            check_true(rank(even2) > 0, "rank(son2_fusion_ring(6)) was not positive")
            check_true(rank(even4) > 0, "rank(son2_fusion_ring(8)) was not positive")

            check_equal(size(multiplication_table(odd)), (rank(odd), rank(odd), rank(odd)),
                "multiplication table size for son2_fusion_ring(5) did not match its rank")
            check_equal(size(multiplication_table(even2)), (rank(even2), rank(even2), rank(even2)),
                "multiplication table size for son2_fusion_ring(6) did not match its rank")
            check_equal(size(multiplication_table(even4)), (rank(even4), rank(even4), rank(even4)),
                "multiplication table size for son2_fusion_ring(8) did not match its rank")
        end

        @testset "basic properties" begin
            check_true(is_commutative(odd), "son2_fusion_ring(5) was not commutative")
            check_true(is_commutative(even2), "son2_fusion_ring(6) was not commutative")
            check_true(is_commutative(even4), "son2_fusion_ring(8) was not commutative")
        end

        @testset "metaplectic alias" begin
            check_equal(rank(meta), rank(odd),
                "metaplectic_fusion_ring(5) did not have the same rank as son2_fusion_ring(5)")
            check_equal(multiplication_table(meta), multiplication_table(odd),
                "metaplectic_fusion_ring(5) did not match son2_fusion_ring(5)")
        end

        @testset "invalid input" begin
            check_throws(
                () -> son2_fusion_ring(3),
                "son2_fusion_ring accepted m=3 even though it should require m ≥ 4"
            )
        end
    end

    #fusion_ring_from_group:

     @testset "fusion_ring_from_group" begin
        z2_tab = [
            1 2
            2 1
        ]

        z3_tab = [
            1 2 3
            2 3 1
            3 1 2
        ]

        bad_not_group = [
            1 1
            2 2
        ]

        @testset "valid group tables" begin
            r2 = fusion_ring_from_group(z2_tab)
            r3 = fusion_ring_from_group(z3_tab)

            check_equal(rank(r2), 2,
                "fusion_ring_from_group on the Z2 Cayley table did not produce rank 2")
            check_equal(rank(r3), 3,
                "fusion_ring_from_group on the Z3 Cayley table did not produce rank 3")

            check_equal(size(multiplication_table(r2)), (2,2,2),
                "fusion_ring_from_group on the Z2 Cayley table had wrong multiplication-table size")
            check_equal(size(multiplication_table(r3)), (3,3,3),
                "fusion_ring_from_group on the Z3 Cayley table had wrong multiplication-table size")

            check_true(is_group_ring(r2),
                "fusion_ring_from_group on the Z2 Cayley table was not detected as a group ring")
            check_true(is_group_ring(r3),
                "fusion_ring_from_group on the Z3 Cayley table was not detected as a group ring")
            check_true(is_commutative(r2),
                "fusion_ring_from_group on the Z2 Cayley table was not commutative")
            check_true(is_commutative(r3),
                "fusion_ring_from_group on the Z3 Cayley table was not commutative")
        end

        @testset "matches zn_fusion_ring on cyclic examples" begin
            r2 = fusion_ring_from_group(z2_tab)
            r3 = fusion_ring_from_group(z3_tab)

            check_equal(multiplication_table(r2), multiplication_table(zn_fusion_ring(2)),
                "fusion_ring_from_group(Z2 table) did not match zn_fusion_ring(2)")
            check_equal(multiplication_table(r3), multiplication_table(zn_fusion_ring(3)),
                "fusion_ring_from_group(Z3 table) did not match zn_fusion_ring(3)")
        end

        @testset "invalid group table is rejected" begin
            check_throws(
                () -> fusion_ring_from_group(bad_not_group),
                "fusion_ring_from_group accepted a table that was not a valid Cayley table"
            )
        end
    end

    #internal group-table helpers:

        @testset "_is_group_table / _inverse_vector" begin
        z2_tab = [
            1 2
            2 1
        ]

        z3_tab = [
            1 2 3
            2 3 1
            3 1 2
        ]

        bad_tab = [
            1 2
            1 2
        ]

        check_true(FusionRings._is_group_table(z2_tab),
            "_is_group_table did not recognize the Z2 Cayley table")
        check_true(FusionRings._is_group_table(z3_tab),
            "_is_group_table did not recognize the Z3 Cayley table")
        check_false(FusionRings._is_group_table(bad_tab),
            "_is_group_table incorrectly accepted an invalid table")

        check_equal(FusionRings._inverse_vector(z2_tab), [1,2],
            "_inverse_vector(Z2 table) was not [1,2]")
        check_equal(FusionRings._inverse_vector(z3_tab), [1,3,2],
            "_inverse_vector(Z3 table) was not [1,3,2]")
    end

    #FusionRing TY

       @testset "FusionRingTY" begin
        z2_tab = [
            1 2
            2 1
        ]

        z3_tab = [
            1 2 3
            2 3 1
            3 1 2
        ]

        bad_tab = [
            1 1
            2 2
        ]

        @testset "basic construction" begin
            r2 = FusionRingTY(z2_tab)
            r3 = FusionRingTY(z3_tab)

            check_equal(rank(r2), 3,
                "FusionRingTY on the Z2 Cayley table did not produce rank 3")
            check_equal(rank(r3), 4,
                "FusionRingTY on the Z3 Cayley table did not produce rank 4")

            check_equal(labels(r2), ["1", "2", "m"],
                "FusionRingTY(Z2) labels were incorrect")
            check_equal(labels(r3), ["1", "2", "3", "m"],
                "FusionRingTY(Z3) labels were incorrect")

            check_equal(size(multiplication_table(r2)), (3,3,3),
                "FusionRingTY(Z2) multiplication table had wrong size")
            check_equal(size(multiplication_table(r3)), (4,4,4),
                "FusionRingTY(Z3) multiplication table had wrong size")
        end

        @testset "basic fusion rules" begin
            r2 = FusionRingTY(z2_tab)

            # last object is m
            m = rank(r2)

            # group element × m = m
            check_equal(fusion_product(r2, 1, m), Dict(m => 1),
                "vacuum × m was not m in FusionRingTY(Z2)")
            check_equal(fusion_product(r2, 2, m), Dict(m => 1),
                "nontrivial group element × m was not m in FusionRingTY(Z2)")

            # m × m = sum of group elements
            check_equal(fusion_product(r2, m, m), Dict(1 => 1, 2 => 1),
                "m × m in FusionRingTY(Z2) was not the sum of all group elements")
        end

        @testset "invalid input" begin
            check_throws(
                () -> FusionRingTY(bad_tab),
                "FusionRingTY accepted an invalid group table"
            )
        end
    end

    #FusiongRingHI:

    
    @testset "FusionRingHI" begin
        z2_tab = [
            1 2
            2 1
        ]

        z3_tab = [
            1 2 3
            2 3 1
            3 1 2
        ]

        nonsymmetric_group_tab = z3_tab  # valid group table, but not symmetric as a matrix
        bad_tab = [
            1 1
            2 2
        ]

        @testset "basic construction on symmetric group table" begin
            r = FusionRingHI(z2_tab)

            check_equal(rank(r), 4,
                "FusionRingHI on the Z2 Cayley table did not produce rank 4")
            check_equal(size(multiplication_table(r)), (4,4,4),
                "FusionRingHI(Z2) multiplication table had wrong size")
            check_equal(labels(r), ["1", "2", "ρ_1", "ρ_2"],
                "FusionRingHI(Z2) labels were incorrect")
        end

        @testset "basic fusion behaviour on Z2 input" begin
            r = FusionRingHI(z2_tab)

            # First two are group sector, last two are rho sector
            check_equal(fusion_product(r, 1, 1), Dict(1 => 1),
                "vacuum × vacuum was not vacuum in FusionRingHI(Z2)")
            check_equal(fusion_product(r, 1, 3), Dict(3 => 1),
                "vacuum × ρ_1 was not ρ_1 in FusionRingHI(Z2)")
        end

        @testset "rejects invalid or unsupported tables" begin
            check_throws(
                () -> FusionRingHI(bad_tab),
                "FusionRingHI accepted a table that was not a valid group table"
            )

            check_throws(
                () -> FusionRingHI(nonsymmetric_group_tab),
                "FusionRingHI accepted a valid group table that was not symmetric as a matrix"
            )
        end
    end

end











