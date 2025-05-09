include("lib/ADMM.jl")
using .ADMM

include("lib/RosenBrock.jl")
import .RosenBrock
include("lib/Booth.jl")
import .Booth

function test1()
    problem = begin
        f(x) = sum(x .^ 2)
        g(z) = sum(abs.(z))
        A = rand(5, 3)
        B = rand(5, 2)
        c = rand(5)
        Problem(f, g, A, B, c)
    end
    println(problem)
    @show admm(problem, zeros(3), zeros(2))
end

function test_rosenbrock()
    problem = begin
        f = RosenBrock.f
        g(z) = sum(abs.(z))
        A = zeros(2, 2)
        B = rand(2, 2)
        c = rand(2)
        Problem(f, g, A, B, c)
    end
    println(problem)
    @show admm(problem, zeros(3), zeros(2))
end
function test_booth()
    problem = begin
        f = Booth.f
        g(z) = sum(abs.(z))
        A = zeros(2, 2)
        B = rand(2, 2)
        c = rand(2)
        Problem(f, g, A, B, c)
    end
    println(problem)
    @show admm(problem, zeros(3), zeros(2))
end
test_rosenbrock()
test_booth()