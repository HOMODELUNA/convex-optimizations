using LinearAlgebra

@enum CondRes Ok = 0 TooBig = 1 TooSmall = -1

"""
计算搜索步长
## 输入
- f 函数
- x,y 当前点的横纵坐标
- d 方向
- current_step 当前步长
- grad_norm 梯度向量的长度
## 返回值
- step 最终得到的步长
- x1,y1 目标点的横纵坐标
"""
function findstep(f, x, y, d, current_step, grad_norm; r=0.001, s=0.618, max_iter=1000)
    a0 = 0 # 初始化搜索范围
    b0 = Inf

    x1 = x + current_step * d
    y1 = f(x1)
    for k = 1:max_iter
        if y1 > y - r * grad_norm * current_step
            b0 = a0
            current_step = a0 + (current_step - a0) * s
        elseif y1 < y - (1 - r) * grad_norm * a0
            b0 = current_step
            current_step = current_step / s
        else
            return current_step, x1, y1
        end
        x1 = x + current_step * d
        y1 = f(x1)
    end
    return current_step, x1, y1
end


# 程序的来源是 https://zhuanlan.zhihu.com/p/44770180
# 解释来自 https://blog.csdn.net/2301_76165902/article/details/142929468
"""
## 输入
- f 一元函数
- g 为f的梯度函数
- x0 初始点
- d 方向
- initial_step_length 初始步长
## 返回值
- α 使 f(x + α)最小
"""
function goldstein(f, g, x0; max_iter=1000, eps=1e-5, initial_step_length=0.001)
    x = x0
    a = initial_step_length
    for k = 1:max_iter
        y = f(x)
        grad = g(x)
        gl = norm(grad)
        if gl^0.5 < eps
            break
        end
        d = -grad / gl
        a, x1, y1 = findstep(f, x, y, d, a, gl)
        if norm(x1 - x, 1) < eps
            break
        end
        x = x1
    end
    return x
end

"""
- a 标量,表示下降步长
- d 矢量, 表示下降方向
- c 控制常数
"""
function satisfy_armijo(f, g, x, d, a; c=0.2)
    f(x + a * d) <= f(x) + c * a * dot(g(x), d)
end

# https://zhuanlan.zhihu.com/p/651901246
"""
# 线搜索回退法
## 输入
- a0 一个足够大的初始步长
"""
function armijo_linear_search_with_retreat(f, g, x, d, a0=100; gamma=0.75, c=0.2)
    a = a0
    while !satisfy_armijo(f, g, x, d, a; c=c)
        a = gamma * a
    end
    return a
end

function satisfy_goldstein(f, g, x, d, a; c_up=0.2, c_down=0.4)
    y0 = f(x)
    x1 = x + a * d
    y1 = f(x1)
    desc = dot(g(x), d)
    if y1 >= y0 + c_up * a * desc
        return TooBig
    elseif y1 >= y0 + c_down * a * desc
        return Ok
    else
        return TooSmall
    end
end

"""
# 线搜索回退法
## 输入
- a0 一个足够大的初始步长
"""
function goldstein_linear_search_with_retreat(f, g, x, d, a0=100; gamma_shrink=0.6, gamma_grow=1.3, c_up=0.2, c_down=0.4)
    a = a0
    while true
        res = satisfy_goldstein(f, g, x, d, a; c_up=c_up, c_down=c_down)
        if res == TooBig
            a = gamma_shrink * a
        elseif res == TooSmall
            a = gamma_grow * a
        else
            return a
        end
    end
    return a
end

# https://www.cnblogs.com/yifdu25/p/8093725.html
function satisfy_powell(f, g, x, d, a; c=0.2, sigma=0.5)
    g0 = g(x)
    x1 = x + a * d
    cond1 = f(x1) <= f(x) + c * a * dot(g0, d)
    return cond1 && dot(g(x1), d) >= sigma * dot(g0, d)
end

"""
# 线搜索回退法
## 输入
- a0 一个足够大的初始步长
"""
function powell_linear_search_with_retreat(f, g, x, d, a0=100; gamma=0.8, c=0.2)
    a = a0
    while !satisfy_powell(f, g, x, d, a; c=c)
        a = gamma * a
    end
    return a
end

function test()
    f(x) = 100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    function g(x)
        d1 = -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1])
        d2 = 200 * (x[2] - x[1]^2)
        return [d1, d2]
    end
    x_0 = [-1, 1]
    d = [1, 1]

    f_a(a) = f(x_0 + d * a)

    @show goldstein(f, g, x_0)
    begin
        method = "armijo"
        a = armijo_linear_search_with_retreat(f, g, x_0, d)
        res = f_a(a)
        println("$method: a=$a, f(x + d * a)=$res")
    end
    begin
        method = "goldstein"
        a = goldstein_linear_search_with_retreat(f, g, x_0, d)
        res = f_a(a)
        println("$method: a=$a, f(x + d * a)=$res")
    end
    begin
        method = "powell"
        a = powell_linear_search_with_retreat(f, g, x_0, d)
        res = f_a(a)
        println("$method: a=$a, f(x + d * a)=$res")
    end
end

test()