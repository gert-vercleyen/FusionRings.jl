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

    