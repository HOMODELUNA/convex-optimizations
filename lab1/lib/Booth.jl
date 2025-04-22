module Booth
export f
"""
最小值在(1,3), 此时 y=0
"""
f(x) = (x[1] + 2 * x[2] - 7)^2 + (2 * x[1] + x[2] - 5)^2

function g(x)
    dx = 2 * (x[1]+2*x[2]-7) + 4*(2*x[1] + x[2]-5)
    dy = 4* (x[1] + 2*x[2]-7) + 2*(2*x[1] + x[2] - 5)
    return [dx,dy]
end

end