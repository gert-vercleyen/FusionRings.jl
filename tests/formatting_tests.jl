@testset "Formatting and printing" begin

    #Integer Transformation helpers

     @testset "transform_integer / bold_integer / subscript_integer / superscript_integer" begin
        check_equal(
            transform_integer(0, bold_digits_dict),
            "𝟎",
            "transform_integer(0, bold_digits_dict) did not return bold zero"
        )

        check_equal(
            transform_integer(123, bold_digits_dict),
            "𝟏𝟐𝟑",
            "transform_integer(123, bold_digits_dict) did not return the expected bold digits"
        )

        check_equal(
            bold_integer(0),
            "𝟎",
            "bold_integer(0) did not return \"𝟎\""
        )

        check_equal(
            bold_integer(907),
            "𝟗𝟎𝟕",
            "bold_integer(907) did not return the expected bold representation"
        )

        check_equal(
            subscript_integer(0),
            "₀",
            "subscript_integer(0) did not return \"₀\""
        )

        check_equal(
            subscript_integer(314),
            "₃₁₄",
            "subscript_integer(314) did not return the expected subscript representation"
        )

        check_equal(
            superscript_integer(0),
            "⁰",
            "superscript_integer(0) did not return \"⁰\""
        )

        check_equal(
            superscript_integer(256),
            "²⁵⁶",
            "superscript_integer(256) did not return the expected superscript representation"
        )
    end

    #Basic string-formatting helpers

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

    @testset "fix_fractions / fix_mult / fix_spaces / fix_powers / fix_cyclo / fix_poly_string" begin
        check_equal(
            fix_fractions("3//4"),
            "\\frac{3}{4}",
            "fix_fractions(\"3//4\") did not convert rational syntax to LaTeX fraction syntax"
        )

        check_equal(
            fix_mult("2*x"),
            "2x",
            "fix_mult(\"2*x\") did not remove multiplication symbols"
        )

        check_equal(
            fix_spaces("a b  c"),
            "abc",
            "fix_spaces did not remove spaces"
        )

        check_equal(
            fix_powers("x^12"),
            "x^{12}",
            "fix_powers(\"x^12\") did not brace a multi-digit exponent"
        )

        check_equal(
            fix_powers("x^2"),
            "x^2",
            "fix_powers(\"x^2\") should not have changed a one-digit exponent"
        )

        check_equal(
            fix_cyclo("zeta(5)"),
            "\\zeta_5",
            "fix_cyclo(\"zeta(5)\") did not convert to zeta subscript notation"
        )

        check_equal(
            fix_poly_string(" 3//4 * x^12 "),
            "\\frac{3}{4}x^{12}",
            "fix_poly_string did not compose the string-cleaning transforms correctly"
        )
    end

    #geometric array helpers:

     @testset "is_geometric_array" begin
        check_true(
            is_geometric_array(Int[]),
            "is_geometric_array(Int[]) should return true"
        )

        check_true(
            is_geometric_array([1]),
            "is_geometric_array([1]) should return true"
        )

        check_true(
            is_geometric_array([1, 2, 4, 8]),
            "is_geometric_array([1,2,4,8]) should return true"
        )

        check_true(
            is_geometric_array([1, -1, 1, -1]),
            "is_geometric_array([1,-1,1,-1]) should return true"
        )

        check_false(
            is_geometric_array([1, 2, 5, 10]),
            "is_geometric_array([1,2,5,10]) should return false"
        )

        check_false(
            is_geometric_array([2]),
            "is_geometric_array([2]) should return false because the one-term case only accepts [1]"
        )
    end

    @testset "factor_squares" begin
        check_equal(
            factor_squares(12),
            (2, 3),
            "factor_squares(12) was not (2,3)"
        )

        check_equal(
            factor_squares(18),
            (3, 2),
            "factor_squares(18) was not (3,2)"
        )

        check_equal(
            factor_squares(16),
            (4, 1),
            "factor_squares(16) was not (4,1)"
        )

        check_equal(
            factor_squares(7),
            (1, 7),
            "factor_squares(7) was not (1,7)"
        )
    end

    #Product_String

     @testset "product_string" begin
        z2 = zn_fusion_ring(2)
        z3 = zn_fusion_ring(3)

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
    end









