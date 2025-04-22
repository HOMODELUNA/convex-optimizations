include("lib/GradientMethods.jl")
using .GradientMethods

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

function test(f,g,title)
    printlnG("--------- function: $title ---------")
    begin
        x = gradient_descend(f, g, [0, 0]; step_length=0.01, max_iter=20000)
        method = "GD"
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end

    begin
        x = gradient_descend_without_zigzag(f, g, [0, 0]; step_length=0.01, max_iter=20000)
        method = "GD de-zigzag"
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    begin
        x = gradient_descend_with_momentum(f, g, [0, 0]; step_length=0.01, max_iter=20000)
        method = "GD momentum"
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    begin
        x = conjugate_gradient(f, g, [0, 0])
        method = "CG "
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    begin
        x = conjugate_gradient_with_momentum(f, g, [0, 0])
        method = "CG momentum"
        y = f(x)
        printlnG("method: $method,\tx = $x,\ty = $(y)")
    end
    println("")
end

test(RosenBrock.f,RosenBrock.g,"RosenBrock")

test(Booth.f,Booth.g,"Booth")