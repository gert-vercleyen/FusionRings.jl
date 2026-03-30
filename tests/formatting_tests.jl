@testset "Formatting and printing" begin

    @testset "integer formatting helpers" begin
        check_equal(
            transform_integer(0, bold_digits_dict),
            "𝟎",
            "transform_integer(0, bold_digits_dict) did not return \"𝟎\""
        )

        check_equal(
            transform_integer(123, bold_digits_dict),
            "𝟏𝟐𝟑",
            "transform_integer(123, bold_digits_dict) did not return the expected bold digits"
        )

        check_equal(
            bold_integer(907),
            "𝟗𝟎𝟕",
            "bold_integer(907) did not return the expected bold-digit string"
        )

        check_equal(
            subscript_integer(314),
            "₃₁₄",
            "subscript_integer(314) did not return the expected subscript string"
        )

        check_equal(
            superscript_integer(256),
            "²⁵⁶",
            "superscript_integer(256) did not return the expected superscript string"
        )
    end

    @testset "element_to_string" begin
        check_equal(
            element_to_string(0, "x"),
            "",
            "element_to_string(0, \"x\") did not return the empty string"
        )

        check_equal(
            element_to_string(1, "x"),
            "x",
            "element_to_string(1, \"x\") did not return \"x\""
        )

        check_equal(
            element_to_string(2, "x"),
            "2 x",
            "element_to_string(2, \"x\") did not return \"2 x\""
        )

        check_equal(
            element_to_string(7, "abc"),
            "7 abc",
            "element_to_string(7, \"abc\") did not return \"7 abc\""
        )
    end

    @testset "polynomial/string cleanup helpers" begin
        check_equal(
            fix_fractions("3//4"),
            "\\frac{3}{4}",
            "fix_fractions(\"3//4\") did not convert to LaTeX fraction syntax"
        )

        check_equal(
            fix_mult("2*x"),
            "2x",
            "fix_mult(\"2*x\") did not remove *"
        )

        check_equal(
            fix_spaces("a b  c"),
            "abc",
            "fix_spaces(\"a b  c\") did not remove spaces"
        )

        check_equal(
            fix_powers("x^12"),
            "x^{12}",
            "fix_powers(\"x^12\") did not brace a multi-digit exponent"
        )

        check_equal(
            fix_powers("x^2"),
            "x^2",
            "fix_powers(\"x^2\") should not change a one-digit exponent"
        )

        check_equal(
            fix_cyclo("zeta(5)"),
            "\\zeta_5",
            "fix_cyclo(\"zeta(5)\") did not convert to zeta subscript notation"
        )

        check_equal(
            fix_poly_string(" 3//4 * x^12 "),
            "\\frac{3}{4}x^{12}",
            "fix_poly_string did not compose the cleanup helpers correctly"
        )
    end

    @testset "factor_squares / is_geometric_array" begin
        check_equal(
            factor_squares(12),
            (2, 3),
            "factor_squares(12) was not (2, 3)"
        )

        check_equal(
            factor_squares(18),
            (3, 2),
            "factor_squares(18) was not (3, 2)"
        )

        check_equal(
            factor_squares(16),
            (4, 1),
            "factor_squares(16) was not (4, 1)"
        )

        check_equal(
            factor_squares(7),
            (1, 7),
            "factor_squares(7) was not (1, 7)"
        )

        check_true(
            is_geometric_array(Int[]),
            "is_geometric_array(Int[]) should be true"
        )

        check_true(
            is_geometric_array([1]),
            "is_geometric_array([1]) should be true"
        )

        check_true(
            is_geometric_array([1, 2, 4, 8]),
            "is_geometric_array([1,2,4,8]) should be true"
        )

        check_true(
            is_geometric_array([1, -1, 1, -1]),
            "is_geometric_array([1,-1,1,-1]) should be true"
        )

        check_false(
            is_geometric_array([1, 2, 5, 10]),
            "is_geometric_array([1,2,5,10]) should be false"
        )
    end

    @testset "product_string" begin
        z2 = zn_fusion_ring(2)
        z3 = zn_fusion_ring(3)
        su2_2 = su2k_fusion_ring(2)

        check_equal(
            product_string(z2, 1, 1),
            "0 × 0 = 0",
            "product_string(Z2, 1, 1) was not \"0 × 0 = 0\""
        )

        check_equal(
            product_string(z2, 2, 2),
            "1 × 1 = 0",
            "product_string(Z2, 2, 2) was not \"1 × 1 = 0\""
        )

        check_equal(
            product_string(z3, 2, 2),
            "1 × 1 = 2",
            "product_string(Z3, 2, 2) was not \"1 × 1 = 2\""
        )

        check_equal(
            product_string(z3, 2, 3),
            "1 × 2 = 0",
            "product_string(Z3, 2, 3) was not \"1 × 2 = 0\""
        )

        check_equal(
            product_string(su2_2, 2, 2),
            "1 × 1 = 0 ⊕ 2",
            "product_string(SU(2)_2, 2, 2) did not show both summands with ⊕"
        )
    end

    @testset "print_multiplication_table / pmt" begin
        z2 = zn_fusion_ring(2)

        out1 = sprint(io -> redirect_stdout(io) do
            print_multiplication_table(z2)
        end)

        check_true(
            occursin("× │ 0 │ 1", out1),
            "print_multiplication_table(Z2) output did not contain the expected header"
        )

        check_true(
            occursin("0 │ 0 │ 1", out1),
            "print_multiplication_table(Z2) output did not contain the expected first row"
        )

        check_true(
            occursin("1 │ 1 │ 0", out1),
            "print_multiplication_table(Z2) output did not contain the expected second row"
        )

        out2 = sprint(io -> redirect_stdout(io) do
            print_multiplication_table(z2; include_zeros=true)
        end)

        check_true(
            occursin("0 0 + 1", out2) || occursin("0 0 + 1 1", out2) || occursin("0 0", out2),
            "print_multiplication_table(Z2; include_zeros=true) did not appear to include zero-multiplicity terms"
        )

        out_pmt = sprint(io -> redirect_stdout(io) do
            pmt(z2)
        end)

        check_equal(
            out_pmt,
            out1,
            "pmt(z2) did not produce the same output as print_multiplication_table(z2)"
        )
    end

    @testset "Base.show(::FusionRing)" begin
        named = zn_fusion_ring(2)
        shown_named = sprint(show, named)

        check_equal(
            shown_named,
            "FR(Z_2)",
            "show(zn_fusion_ring(2)) did not use the ring name branch"
        )

        mt = zeros(Int, 2, 2, 2)
        mt[1,1,1] = 1
        mt[1,2,2] = 1
        mt[2,1,2] = 1
        mt[2,2,1] = 1

        unnamed = fusion_ring(mt; labels=["a","b"])
        shown_unnamed = sprint(show, unnamed)

        check_true(
            startswith(shown_unnamed, "FR("),
            "show on an unnamed fusion ring did not begin with \"FR(\""
        )

        check_true(
            endswith(shown_unnamed, ")"),
            "show on an unnamed fusion ring did not end with \")\""
        )
    end

    @testset "default labels use bold_integer" begin
        mt = zeros(Int, 2, 2, 2)
        mt[1,1,1] = 1
        mt[1,2,2] = 1
        mt[2,1,2] = 1
        mt[2,2,1] = 1

        r = fusion_ring(mt)

        check_equal(
            labels(r)[1],
            bold_integer(1),
            "fusion_ring(mt) did not use bold_integer(1) as the first default label"
        )

        check_equal(
            labels(r)[2],
            bold_integer(2),
            "fusion_ring(mt) did not use bold_integer(2) as the second default label"
        )
    end

end