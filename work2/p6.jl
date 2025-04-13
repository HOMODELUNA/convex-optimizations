using LinearAlgebra
# https://zhuanlan.zhihu.com/p/306635632


"""
# DFP 拟牛顿法
## 输入
- f: 一元函数
- gradient: f 的梯度
- x_init : 初始值
- eps: 计算精度
## 返回值
一个最小值点
"""
function broyden_flether_goldfarb_shanno(f,gradient,x_init; eps=0.03)
    H_0 = (l = length(x_init); Matrix(I,l,l))
    g_0 = gradient(x_init)
    if norm(g_0) < eps 
        return x_init
    end
    function iter(k,x,H,g)
        if (k <=0 ) 
            return x
        end
        d_k = - H * g
        lambda_k = calculate_proper_steplength(f,x,d_k)
        s_k = lambda_k * d_k
        x_next = x + s_k
        
        g_next = gradient(x_next)
        if norm(g_next) < eps
            return x_next
        end
        y = g_next - g
        H_next = H + (s_k * transpose(s_k)) / dot(s_k,y) - (H * y * transpose(y) * H) / (transpose(y) * H * y)
        return iter(k-1,x_next,H_next,g_next)
    end
    return iter(100,x_init,H_0,g_0)
end

# BFGS optimization algorithm
function bfgs(f,gradient,x0; eps=1e-6, max_iter=1000)
    n = length(x0)
    H = I(n)  # Initial Hessian approximation
    x = x0

    for iter in 1:max_iter
        g = gradient(x)
        if norm(g) < eps
            println("Converged after $iter iterations.")
            return x
        end
        
        # Search direction
        p = -H * g
        
        # Line search (simple step)
        alpha = 0.1  # Step size
        x_new = x + alpha * p
        
        # Compute new gradient
        g_new = gradient(x_new)
        
        # Compute the difference in position and gradient
        s = x_new - x
        y = g_new - g

        # Update Hessian approximation using BFGS formula
        rho = 1 / (y' * s)
        H = (I(n) - rho * s * y') * H * (I(n) - rho * y * s') + rho * s * s'
        
        # Update x
        x = x_new
    end

    println("Maximum iterations reached.")
    return x
end

function test_bfgs()
f(x) = x[1]^2 + 4 * x[2]^2 - 4*x[1]-8*x[2]
function g(x) 
     [2 * x[1] - 4,8  * x[2] - 8]
end
H_0 = I
x_0 = [0,0]
x_min =  bfgs(f,g,x_0)
y = f(x_min)
println("BFGS x_min = $x_min, y=$y")
end

test_bfgs()