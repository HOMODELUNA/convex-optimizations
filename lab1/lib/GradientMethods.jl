module GradientMethods
include("LinearSearch.jl")
export gradient_descend
export gradient_descend_without_zigzag
export gradient_descend_with_momentum
export conjugate_gradient
export conjugate_gradient_with_momentum
using LinearAlgebra
import .LinearSearch
"""
梯度下降法
## 参数
- f(x) : 一阶函数
- g(x) : f的梯度函数
- x0 : 初始值
- step_length :  步长
- eps : 误差限
- max_iter : 最大迭代次数
## 返回值
- 最小值点
"""
function gradient_descend(f, g, x0; step_length=0.05, eps=0.001, max_iter=10000)
    x = x0
    for k = 1:max_iter
        grad = g(x)
        n_grad = norm(grad)
        if n_grad < eps
            println("-- GD return at $k iter")
            return x
        end
        step = (grad / n_grad) * step_length
        x1 = x - step
        x = x1
    end
    println("-- GD max iter $max_iter retrieved")
    return x
end

function gradient_descend_without_zigzag(f, g, x0; step_length=0.05, eps=0.001, max_iter=10000)
    x = x0
    grad_prev = g(x)
    for k = 1:max_iter
        grad = g(x)
        if dot(grad, grad_prev) < 0 # 拐弯幅度过大
            grad -= grad_prev * dot(grad_prev, grad) / dot(grad_prev, grad_prev)
        end
        n_grad = norm(grad)
        if n_grad < eps
            return x
        end
        step = (grad / n_grad) * step_length
        x1 = x - step
        x = x1
    end
    println("-- GD de-zigzag max iter $max_iter retrieved")
    return x
end

function gradient_descend_with_momentum(f, g, x0; step_length=0.05, eps=0.001, max_iter=10000)
    x = x0
    momentum = zero(x)
    for k = 1:max_iter
        grad = g(x)
        n_grad = norm(grad)
        if n_grad < eps
            println("-- GD momentum return at $k iter")
            return x
        end
        step = (grad / n_grad) * step_length
        momentum = momentum * 0.9 + (-step * 0.1)
        x1 = x - step + momentum
        x = x1
    end
    println("-- GD momentum max iter $max_iter retrieved")
    return x
end

function get_min_alpha(f, g, x, d; eps=0.00001)
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

"""
通过Crowder_wolfe公式计算beta值
## 返回值: beta
"""
function crowder_wolfe(g1, d0, G)
    return (transpose(g1) * G * d0) / (transpose(d0) * G * d0)
end

"""
共轭方向法
"""
function conjugate_gradient(f, g, x0; eps=0.01, max_iter=10000)
    G = I(2)
    g0 = g(x0)
    # 循环变量
    d = g0
    x = x0
    for k = 1:max_iter
        if norm(d) < eps
            println("-- CG return at $k iter")
            return x
        end
        alpha = get_min_alpha(f, g, x, -d)
        @assert !isnan(alpha)
        x1 = x - alpha * d

        g1 = g(x1)

        beta = crowder_wolfe(g1, d, G)
        d1 = -g1 + beta * d
        # if k <= 20
        #     println("iter $k, alpha = $alpha,d  = $d => $d1, g1 = $g1, x = $x => $x1")
        # end
        # if k % 100 == 0
        #     println("iter $k, alpha = $alpha,d  = $d => $d1, g1 = $g1, x = $x => $x1")
        # end

        # 下个周期
        x = x1
        d = d1
    end
    println("-- CG max iter $max_iter retrieved")
    return x
end
function conjugate_gradient_with_momentum(f, g, x0; eps=0.01, max_iter=10000, momentum_factor=0.1)
    G = I(2)
    g0 = g(x0)
    momentum = zero(g0)
    # 循环变量
    d = g0
    x = x0
    for k = 1:max_iter
        if norm(d) < eps
            println("-- CG Momentum return at $k iter")
            return x
        end
        alpha = get_min_alpha(f, g, x, -d)
        @assert !isnan(alpha)
        x1 = x - alpha * d

        g1 = g(x1)
        if norm(g1) < eps
            println("-- CG Momentum return at $k iter")
            return x
        end
        beta = crowder_wolfe(g1, d, G)
        d1 = -g1 + beta * d + 0.5 * momentum
        momentum = (1-momentum_factor)* momentum + momentum_factor * d1
        # if k <= 20
        #     println("iter $k, alpha = $alpha,d  = $d => $d1, g1 = $g1, x = $x => $x1")
        # end
        # if k % 100 == 0
        #     println("iter $k, alpha = $alpha,d  = $d => $d1, g1 = $g1, x = $x => $x1")
        # end

        # 下个周期
        x = x1
        d = d1
    end
    println("CG momentum max iter $max_iter retrieved")
    return x
end
end