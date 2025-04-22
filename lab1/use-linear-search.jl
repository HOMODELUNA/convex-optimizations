include("lib/linear_search.jl")

using .LinearSearch


test_1(x) = 2x^2 - x - 1
test_2(x) = 3x^2 - 21.6 * x - 1

@show golden_section(test_1, -1, 1, 0.06)
@show golden_section(test_2, 0, 25, 0.08)

@show fibonacci_subsequence(test_1, -1, 1, 0.06)
@show fibonacci_subsequence(test_2, 0, 25, 0.08)

@show bisection(test_1, -1, 1, 0.06)
@show bisection(test_2, 0, 25, 0.08)

@show shubert_piyavskii(test_1, -1, 1, 0.06)
@show shubert_piyavskii(test_2, 0, 25, 0.08)