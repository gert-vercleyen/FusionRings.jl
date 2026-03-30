@testset "Formatting and printing" begin

    # ------------------------------------------------------------
    # transform_integer / bold_integer / subscript_integer / superscript_integer
    # ------------------------------------------------------------

    @testset "integer formatting helpers" begin
        @testset "single digits" begin
            check_equal(transform_integer(0, bold_digits_dict), "𝟎",
                "transform_integer(0, bold_digits_dict) did not return bold zero")
            check_equal(transform_integer(5, bold_digits_dict), "𝟓",
                "transform_integer(5, bold_digits_dict) did not return bold five")

            check_equal(bold_integer(0), "𝟎",
                "bold_integer(0) did not return \"𝟎\"")
            check_equal(subscript_integer(0), "₀",
                "subscript_integer(0) did not return \"₀\"")
            check_equal(superscript_integer(0), "⁰",
                "superscript_integer(0) did not return \"⁰\"")
        end

        @testset "multiple digits" begin
            check_equal(transform_integer(123, bold_digits_dict), "𝟏𝟐𝟑",
                "transform_integer(123, bold_digits_dict) did not return the expected bold string")
            check_equal(bold_integer(907), "𝟗𝟎𝟕",
                "bold_integer(907) did not return the expected bold string")
            check_equal(subscript_integer(314), "₃₁₄",
                "subscript_integer(314) did not return the expected subscript string")
            check_equal(superscript_integer(256), "²⁵⁶",
                "superscript_integer(256) did not return the expected superscript string")
        end

        @testset "repeated digits / internal zeros" begin
            check_equal(bold_integer(1001), "𝟏𝟎𝟎𝟏",
                "bold_integer(1001) did not preserve repeated digits and zeros")
            check_equal(subscript_integer(1001), "₁₀₀₁",
                "subscript_integer(1001) did not preserve repeated digits and zeros")
            check_equal(superscript_integer(1001), "¹⁰⁰¹",
                "superscript_integer(1001) did not preserve repeated digits and zeros")
        end

        @testset "string return type" begin
            check_true(bold_integer(12) isa AbstractString,
                "bold_integer(12) did not return a string")
            check_true(subscript_integer(12) isa AbstractString,
                "subscript_integer(12) did not return a string")
            check_true(superscript_integer(12) isa AbstractString,
                "superscript_integer(12) did not return a string")
        end
    end

    # ------------------------------------------------------------
    # element_to_string
    # ------------------------------------------------------------

    @testset "element_to_string" begin
        @testset "zero multiplicity" begin
            check_equal(element_to_string(0, "x"), "",
                "element_to_string(0, \"x\") did not return the empty string")
            check_equal(element_to_string(0, "abc"), "",
                "element_to_string(0, \"abc\") did not return the empty string")
        end

        @testset "unit multiplicity" begin
            check_equal(element_to_string(1, "x"), "x",
                "element_to_string(1, \"x\") did not return \"x\"")
            check_equal(element_to_string(1, "abc"), "abc",
                "element_to_string(1, \"abc\") did not return \"abc\"")
        end

        @testset "higher multiplicity" begin
            check_equal(element_to_string(2, "x"), "2 x",
                "element_to_string(2, \"x\") did not return \"2 x\"")
            check_equal(element_to_string(7, "abc"), "7 abc",
                "element_to_string(7, \"abc\") did not return \"7 abc\"")
        end
    end

    # cleanup helpers

    @testset "fix_fractions" begin
        check_equal(fix_fractions("3//4"), "\\frac{3}{4}",
            "fix_fractions(\"3//4\") did not convert a simple rational")

        check_equal(fix_fractions("1//2 + 3//4"), "\\frac{1}{2} + \\frac{3}{4}",
            "fix_fractions did not convert multiple rationals in one string")

        check_equal(fix_fractions("x"), "x",
            "fix_fractions should leave strings without // unchanged")
    end

    @testset "fix_mult" begin
        check_equal(fix_mult("2*x"), "2x",
            "fix_mult(\"2*x\") did not remove *")
        check_equal(fix_mult("a*b*c"), "abc",
            "fix_mult(\"a*b*c\") did not remove all *")
        check_equal(fix_mult("xyz"), "xyz",
            "fix_mult should leave strings without * unchanged")
    end

    @testset "fix_spaces" begin
        check_equal(fix_spaces("a b  c"), "abc",
            "fix_spaces(\"a b  c\") did not remove spaces")
        check_equal(fix_spaces("   x   "), "x",
            "fix_spaces(\"   x   \") did not remove leading/trailing spaces")
        check_equal(fix_spaces("xyz"), "xyz",
            "fix_spaces should leave strings without spaces unchanged")
    end

    @testset "fix_powers" begin
        check_equal(fix_powers("x^12"), "x^{12}",
            "fix_powers(\"x^12\") did not brace a multi-digit exponent")
        check_equal(fix_powers("x^2"), "x^2",
            "fix_powers(\"x^2\") should not change a one-digit exponent")
        check_equal(fix_powers("x^12 + y^34"), "x^{12} + y^{34}",
            "fix_powers did not handle multiple multi-digit exponents")
        check_equal(fix_powers("abc"), "abc",
            "fix_powers should leave strings without exponents unchanged")
    end

    @testset "fix_cyclo" begin
        check_equal(fix_cyclo("zeta(5)"), "\\zeta_5",
            "fix_cyclo(\"zeta(5)\") did not convert to zeta subscript notation")
        check_equal(fix_cyclo("zeta(3) + zeta(7)"), "\\zeta_3 + \\zeta_7",
            "fix_cyclo did not handle multiple zeta terms")
        check_equal(fix_cyclo("x"), "x",
            "fix_cyclo should leave strings without zeta(...) unchanged")
    end

    @testset "fix_poly_string" begin
        check_equal(fix_poly_string(" 3//4 * x^12 "), "\\frac{3}{4}x^{12}",
            "fix_poly_string did not compose cleanup helpers correctly")

        check_equal(
            fix_poly_string(" zeta(5) * x^12 + 1//2 * y^34 "),
            "\\zeta_5x^{12}+\\frac{1}{2}y^{34}",
            "fix_poly_string did not correctly combine cyclotomic, fraction, multiplication, space, and power cleanup"
        )

        check_equal(fix_poly_string("x"), "x",
            "fix_poly_string should leave already-simple strings unchanged")
    end

    # factor_squares

    @testset "factor_squares" begin
        check_equal(factor_squares(1), (1, 1),
            "factor_squares(1) was not (1,1)")
        check_equal(factor_squares(7), (1, 7),
            "factor_squares(7) was not (1,7)")
        check_equal(factor_squares(12), (2, 3),
            "factor_squares(12) was not (2,3)")
        check_equal(factor_squares(18), (3, 2),
            "factor_squares(18) was not (3,2)")
        check_equal(factor_squares(16), (4, 1),
            "factor_squares(16) was not (4,1)")
        check_equal(factor_squares(72), (6, 2),
            "factor_squares(72) was not (6,2)")
    end

    # is_geometric_array

    @testset "is_geometric_array" begin
        @testset "trivial cases" begin
            check_true(is_geometric_array(Int[]),
                "is_geometric_array(Int[]) should be true")
            check_true(is_geometric_array([1]),
                "is_geometric_array([1]) should be true")
            check_false(is_geometric_array([2]),
                "is_geometric_array([2]) should be false")
        end

        @testset "positive ratio" begin
            check_true(is_geometric_array([1,2,4,8]),
                "is_geometric_array([1,2,4,8]) should be true")
            check_true(is_geometric_array([3,6,12,24]),
                "is_geometric_array([3,6,12,24]) should be true")
        end

        @testset "negative ratio" begin
            check_true(is_geometric_array([1,-1,1,-1]),
                "is_geometric_array([1,-1,1,-1]) should be true")
            check_true(is_geometric_array([2,-4,8,-16]),
                "is_geometric_array([2,-4,8,-16]) should be true")
        end

        @testset "non-geometric" begin
            check_false(is_geometric_array([1,2,5,10]),
                "is_geometric_array([1,2,5,10]) should be false")
            check_false(is_geometric_array([1,3,9,28]),
                "is_geometric_array([1,3,9,28]) should be false")
        end
    end

    # product_string

    @testset "product_string" begin
        z2 = zn_fusion_ring(2)
        z3 = zn_fusion_ring(3)
        su2_2 = su2k_fusion_ring(2)

        @testset "pointed rings" begin
            check_equal(product_string(z2, 1, 1), "0 × 0 = 0",
                "product_string(Z2,1,1) was not \"0 × 0 = 0\"")
            check_equal(product_string(z2, 2, 2), "1 × 1 = 0",
                "product_string(Z2,2,2) was not \"1 × 1 = 0\"")
            check_equal(product_string(z3, 2, 2), "1 × 1 = 2",
                "product_string(Z3,2,2) was not \"1 × 1 = 2\"")
            check_equal(product_string(z3, 2, 3), "1 × 2 = 0",
                "product_string(Z3,2,3) was not \"1 × 2 = 0\"")
        end

        @testset "multiple summands" begin
            check_equal(product_string(su2_2, 2, 2), "1 × 1 = 0 ⊕ 2",
                "product_string(SU(2)_2,2,2) did not show both summands correctly")
        end

        @testset "output shape" begin
            s = product_string(z3, 2, 2)
            check_true(occursin(" × ", s),
                "product_string output did not contain the multiplication symbol")
            check_true(occursin(" = ", s),
                "product_string output did not contain the equals sign")
        end
    end

    # print_multiplication_table / pmt

    @testset "print_multiplication_table / pmt" begin
        z2 = zn_fusion_ring(2)
        z3 = zn_fusion_ring(3)

        @testset "smoke tests" begin
            ok1 = true
            try
                sprint(io -> redirect_stdout(io) do
                    print_multiplication_table(z2)
                end)
            catch err
                ok1 = false
                println("FAILED: print_multiplication_table(Z2) threw an error")
                println("  error: ", sprint(showerror, err))
            end
            @test ok1

            ok2 = true
            try
                sprint(io -> redirect_stdout(io) do
                    print_multiplication_table(z3; include_zeros=true)
                end)
            catch err
                ok2 = false
                println("FAILED: print_multiplication_table(Z3; include_zeros=true) threw an error")
                println("  error: ", sprint(showerror, err))
            end
            @test ok2
        end

        @testset "contains labels" begin
            out = sprint(io -> redirect_stdout(io) do
                print_multiplication_table(z2)
            end)

            check_true(occursin("0", out),
                "print_multiplication_table(Z2) output did not contain label 0")
            check_true(occursin("1", out),
                "print_multiplication_table(Z2) output did not contain label 1")
            check_true(occursin("×", out),
                "print_multiplication_table(Z2) output did not contain the multiplication symbol")
        end

        @testset "alias agreement" begin
            out1 = sprint(io -> redirect_stdout(io) do
                print_multiplication_table(z2)
            end)

            out2 = sprint(io -> redirect_stdout(io) do
                pmt(z2)
            end)

            check_equal(out2, out1,
                "pmt(z2) did not match print_multiplication_table(z2)")
        end
    end

    # Base.show(::FusionRing)

    @testset "Base.show(::FusionRing)" begin
        @testset "named ring branch" begin
            z2 = zn_fusion_ring(2)
            shown = sprint(show, z2)

            check_true(startswith(shown, "FR("),
                "show(zn_fusion_ring(2)) did not begin with \"FR(\"")
            check_true(occursin("Z_2", shown),
                "show(zn_fusion_ring(2)) did not contain the ring name")
            check_true(endswith(shown, ")"),
                "show(zn_fusion_ring(2)) did not end with \")\"")
        end

        @testset "unnamed ring fallback branch" begin
            mt = zeros(Int, 2, 2, 2)
            mt[1,1,1] = 1
            mt[1,2,2] = 1
            mt[2,1,2] = 1
            mt[2,2,1] = 1

            r = fusion_ring(mt; labels=["a","b"])
            shown = sprint(show, r)

            check_true(startswith(shown, "FR("),
                "show on an unnamed ring did not begin with \"FR(\"")
            check_true(endswith(shown, ")"),
                "show on an unnamed ring did not end with \")\"")
        end
    end

    # integration: default labels
    # 

    @testset "default labels use bold_integer" begin
        mt = zeros(Int, 2, 2, 2)
        mt[1,1,1] = 1
        mt[1,2,2] = 1
        mt[2,1,2] = 1
        mt[2,2,1] = 1

        r = fusion_ring(mt)

        check_equal(labels(r)[1], bold_integer(1),
            "fusion_ring(mt) did not use bold_integer(1) as the first default label")
        check_equal(labels(r)[2], bold_integer(2),
            "fusion_ring(mt) did not use bold_integer(2) as the second default label")
    end

end