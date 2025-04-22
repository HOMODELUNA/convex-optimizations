module NoGradientMethods
export cyclic_coordinate
export cyclic_coordinate_with_momentum
export hooke_jeeves

include("LinearSearch.jl")
import .LinearSearch
using LinearAlgebra

"""
坐标轮换法
## 输入
- f 一元函数
- x0 初始值
## 返回值
x 最小值点
"""
function cyclic_coordinate(f, x0; eps=0.001, max_iter=10000)
    coordinates = length(x0)
    x = copy(x0)
    for k = 1:max_iter
        x1 = copy(x)
        for c = 1:coordinates
            function f_on_coordinate(d)
                x2 = copy(x1)
                x2[c] += d
                return f(x2)
            end
            res = LinearSearch.golden_section(f_on_coordinate, -50.0, 50.0, 0.0001)
            x1[c] += res.point
        end
        if norm(x1 - x) < eps
            println("-- Cyclic Coordinate return at iter $max_iter")
            return x1
        end
        x = x1
    end
    println("-- Cyclic Coordinate max iter $max_iter exeeced")
    return x
end
"""
坐标轮换法(添加动量因素)
## 输入
- f 一元函数
- x0 初始值
## 返回值
x 最小值点
"""
function cyclic_coordinate_with_momentum(f, x0; eps=0.001, max_iter=10000)
    coordinates = length(x0)
    x = copy(x0)
    momentum = zero(x0)
    for k = 1:max_iter
        x1 = copy(x)
        for c = 1:coordinates
            function f_on_coordinate(d)
                x2 = copy(x1)
                x2[c] += d
                return f(x2)
            end
            res = LinearSearch.golden_section(f_on_coordinate, -50.0, 50.0, 0.0001)
            x1[c] += res.point
        end
        if norm(x1 - x) < eps
            println("-- Cyclic Coordinate Momentum return at iter $max_iter")
            return x1
        end
        momentum = 0.9 * momentum + 0.1 * (x1 - x)
        x = x1 + momentum
    end
    println("-- Cyclic Coordinate Momentum max iter $max_iter exeeced")
    return x
end

# 这个讲得最好
# https://web.stanford.edu/group/sisl/k12/optimization/MO-unit2-pdfs/2.11minimum3D2hooke-jeeves.pdf
"""
- beta : 放大系数
- alpha : 缩小系数
"""
function hooke_jeeves(f, x0; eps=0.001, beta=2.0, step_size=0.1, alpha=0.25, max_iter=10000)
    x = x0
    for k = 1:max_iter
        x_improve = hj_compare(f, x, step_size)
        if x_improve == x # 没有改进方向,收缩步长
            if step_size < eps 
                println("-- Hooke & Jeeves return at iter $k iter")
                return x
            end
            step_size *= alpha
            continue
        end
        # 有改进方向,沿该方向延申
        d = x_improve - x
        x_extend = hj_go_through(f, x_improve, d)
        x = x_extend
    end
    println("-- Hooke & Jeeves Momentum max iter $max_iter exeeced")
    return x
end

function direction(x, c)
    e = zero(x)
    e[c] = 1.0
    return e
end



function hj_compare(f, x0, delta)
    function try_coord(e, x_old)
        y_old = f(x_old)
        x_new = x_old + delta * e
        y_new = f(x_new)
        if y_new < y_old
            return x_new
        end
        x_new = x_old - delta * e
        y_new = f(x_new)
        if y_new < y_old
            return x_new
        end
        return x_old
    end
    coords = length(x0)
    x = x0
    for c = 1:coords
        e = direction(x0, c)
        x = try_coord(e, x)
    end
    return x
end

function hj_go_through(f, x0, d)
    x = x0
    y = f(x)
    while true
        x1 = x + d
        y1 = f(x1)
        if y1 > y
            return x
        end
        y = y1
        x = x1
    end
end
end