include("lib/rosenbrock.jl")
import .RosenBrock

@show RosenBrock.f([0, 1])
@show RosenBrock.g([1, 1])