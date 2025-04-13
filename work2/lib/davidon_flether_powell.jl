using LinearAlgebra

"""
### 返回值

(l,r) 一个区间,相信其中含有最小值
"""
function one_d_min_try(f,x,delta)
    f_0 = f(x)
    h = 0.5
    
    f_l = f(x - h * delta)
    f_r = f(x + h * delta)
    if (f_l > f_0 && f_r > f_0)
        return -h,h
    end
    if (f_l > f_0 && f_0 > f_r)
        # 正向查找
        function iter((d0,f0), (d1,f1))
            d2 = 2*d1
            
            f2 = f(x + d2 * delta)
            if f2 > f1
                return d0,d2
            end
            return iter((d1,f1),(d2,f2))
        end
        return iter((0,f_0),(h,f_r))
    end
    if (f_l < f_0 && f_0 < f_r)
        # 负向查找
        (l, r) = one_d_min_try(f,x,-delta)
        return (-r,-l)
    end
    error("the function may not be convex")
end

include("fibbonacci_subsequence.jl")
function test_one_d_min_try() 
    test1(x) = (x-1)*(x-3)
    @show one_d_min_try(test1,0,1)
    (l,r) = one_d_min_try(test1,0,1)

    @show fibonacci_subsequence(test1,l,r,0.03)
end


function calculate_proper_steplength(f,x,delta, eps=0.03)
    (l,r) = one_d_min_try(f,x,delta)
    g(d) = f( x + d * delta)
    (; point , value) = fibonacci_subsequence(g,l,r,eps)
    return point
end

# https://zhuanlan.zhihu.com/p/32094294
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
function davidon_flether_powell(f,gradient,x_init; eps=0.03)
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