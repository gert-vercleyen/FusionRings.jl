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


