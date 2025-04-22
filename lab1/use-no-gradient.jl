include("lib/NoGradientMethods.jl")
using .NoGradientMethods

include("lib/RosenBrock.jl")
import .RosenBrock

include("lib/Booth.jl")
import .Booth

GREEN = "\033[32m"
FINISH = "\033[0m"
function printlnG(x)
    print(GREEN)
    print(x)
    println(FINISH)
end

function test(f,title)
    printlnG("--------- function: $title ---------")
    begin
        x = cyclic_coordinate(f, [0.0, 0.0]; max_iter=20000)
        method = "Cyclic Coordinate"
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    begin
        x = cyclic_coordinate_with_momentum(f, [0.0, 0.0]; max_iter=20000)
        method = "Cyclic Coordinate Momentum"
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    begin
        x = hooke_jeeves(f, [0.0, 0.0]; eps=0.001,max_iter=20000)
        method = "Hooke & Jeeves"
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    println("")
end

test(RosenBrock.f,"RosenBrock")


function f2(arr)
x = arr[1];
y = arr[2]
return 2*(x-1) ^2 + (y-2)^2
end

test(f2,"f(x) = 2(x[1]-1)^2 + (x[2]-2)^2")
test(Booth.f,"Booth")
