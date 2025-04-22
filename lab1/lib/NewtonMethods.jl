module NewtonMethods
export newton,newton_damping
using LinearAlgebra
include("LinearSearch.jl")
import .LinearSearch
# https://www.cnblogs.com/MarisaMagic/p/17904136.html
"""
牛顿法
## 输入
- f(x) 一元函数
- g(x) f的梯度
- g2(x) 计算f的Hessian 矩阵 
- x0 初始值
"""
function newton(f,g,g2,x0; eps = 0.01,step = 1,max_iter=10000)
    x = x0
    for k = 1:max_iter
        grad = g(x)
        if norm(grad) < eps
            println("-- Newton return at $k iter")
            return x
        end
        x1 = x - step * inv(g2(x)) * grad
        x = x1
    end
    println("-- Newton max iter $max_iter exceeded")
    return x
end

"""
阻尼牛顿法,在牛顿法的基础上,通过线搜索确定每步的步长
## 输入
- f(x) 一元函数
- g(x) f的梯度
- g2(x) 计算f的Hessian 矩阵 
- x0 初始值
"""
function newton_damping(f,g,g2,x0; eps = 0.01,step = 1,max_iter=10000)
    x = x0
    for k = 1:max_iter
        grad = g(x)
        if norm(grad) < eps
            println("-- Newton return at $k iter")
            return x
        end
        d = - (inv(g2(x)) * grad)
        step = get_min_alpha(f,x,d)
        x1 = x + step * d
        x = x1
    end
    println("-- Newton max iter $max_iter exceeded")
    return x
end

function get_min_alpha(f, x, d; eps=0.00001)
    inner_f(a) = f(x + a * d)
    y = f(x)
    res = LinearSearch.golden_section(inner_f, 0, 10, eps)
    if res.value <= y
        return res.point
    else
    end
    res2 = LinearSearch.golden_section(inner_f, -10, 0, eps)
    if res2.value < y
        return res2.point
    end

    res3 = LinearSearch.golden_section(inner_f, res2.point, res.point, eps / 10)
    if res3.value > f(x)
        error("invalid alpha, x0 = $x, alpha=$(res.point), y = $(f(x)), next_y = $(res.value)")
    end
    # println("get min alpha, x0 = $x, alpha=$(res.point), y = $(f(x)), next_y = $(res.value)")
    return res.point
end
end