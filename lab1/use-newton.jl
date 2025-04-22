include("lib/NewtonMethods.jl")
using .NewtonMethods

include("lib/RosenBrock.jl")

GREEN = "\033[32m"
FINISH = "\033[0m"
function printlnG(x)
    print(GREEN)
    print(x)
    println(FINISH)
end

begin
    using .RosenBrock
    begin
        method = "Newton"
        x = newton(f, g, g2, [0, 0])
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    begin
        method = "Newton Damping"
        x = newton_damping(f, g, g2, [0, 0])
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
end