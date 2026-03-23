using Test

function check_true(cond::Bool, msg::AbstractString)
    cond || println("FAILED: ", msg)
    @test cond
end

function check_false(cond::Bool, msg::AbstractString)
    cond && println("FAILED: ", msg)
    @test !cond
end

function check_equal(actual, expected, msg::AbstractString)
    if actual != expected
        println("FAILED: ", msg)
        println("  expected: ", repr(expected))
        println("  got     : ", repr(actual))
    end
    @test actual == expected
end

function check_approx(actual, expected, msg::AbstractString; atol=1e-8, rtol=1e-8)
    ok = isapprox(actual, expected; atol=atol, rtol=rtol)
    if !ok
        println("FAILED: ", msg)
        println("  expected approximately: ", repr(expected))
        println("  got                   : ", repr(actual))
    end
    @test isapprox(actual, expected; atol=atol, rtol=rtol)
end

function check_throws(f::Function, msg::AbstractString)
    threw = false
    try
        f()
    catch
        threw = true
    end
    threw || println("FAILED: ", msg)
    @test threw
end