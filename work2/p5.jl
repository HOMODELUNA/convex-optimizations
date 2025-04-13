include("lib/davidon_flether_powell.jl")

begin
    function test_f(x)
        x1 = x[1]
        x2 = x[2]
        return 10 * x1^2 + x2^2
    end
    function test_g(x)
        x1 = x[1]
        x2 = x[2]
        return [20 * x1 ,  2 * x2]
    end
    x_0 = [0.1,1]
    @show davidon_flether_powell(test_f,test_g,x_0; eps= 0.03)
end