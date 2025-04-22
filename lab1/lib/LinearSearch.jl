module LinearSearch
# https://juliacollections.github.io/DataStructures.jl/latest/
using DataStructures
export Solution, golden_section, bisection, shubert_piyavskii,fibonacci_subsequence
struct Solution
    point
    value
end

λ = (sqrt(5) - 1) / 2

# https://zhuanlan.zhihu.com/p/666835111
"""
黄金分割法
## 参数
- f(x) 一元函数
- l 区间左值
- r 区间右值
- eps 误差限
## 返回值
Solution 结构体,含有
- point 最小值点
- value 最小值
"""
function golden_section(f, l, r, eps)
    if l > r
        l, r = r, l
    end
    if r - l < eps
        p = (l + r) / 2
        return Solution(p, f(p))
    end
    d = r - l
    m_right = l + λ * d
    m_left = r - λ * d
    if f(m_right) > f(m_left)
        return golden_section(f, l, m_right, eps)
    else
        return golden_section(f, m_left, r, eps)
    end
end

# https://blog.csdn.net/tianhai12/article/details/130491491
"""
斐波那契数列法
## 参数
- f(x) 一元函数
- l 区间左值
- r 区间右值
- eps 误差限
## 返回值
Solution 结构体,含有
- point 最小值点
- value 最小值
"""
function fibonacci_subsequence(f, l, r, eps, f_n0=1, f_n1=1, f_n2=2)
    if l > r
        l, r = r, l
    end
    if r - l < eps
        p = (l + r) / 2
        return Solution(p, f(p))
    end
    d = r - l
    m_left = l + f_n0 * d / f_n2
    m_right = l + f_n1 * d / f_n2
    if f(m_left) > f(m_right)
        return fibonacci_subsequence(f, m_left, r, eps, f_n1, f_n2, f_n1 + f_n2)
    else
        return fibonacci_subsequence(f, l, m_right, eps, f_n1, f_n2, f_n1 + f_n2)
    end
end

"""
二分法
## 参数
- f(x) 一元函数
- l 区间左值
- r 区间右值
- eps 误差限
## 返回值
Solution 结构体,含有
- point 最小值点
- value 最小值
"""
# https://blog.csdn.net/qq_36089856/article/details/103915285
function bisection(f, l, r, eps)
    if l > r
        l, r = r, l
    end
    if r - l < eps
        p = (l + r) / 2
        return Solution(p, f(p))
    end
    mid = (l + r) / 2
    eta = min(eps, (r - l) / 4)
    m_l = mid - eta
    m_r = mid + eta
    if f(m_l) > f(m_r)
        return bisection(f, m_l, r, eps)
    else
        return bisection(f, l, m_r, eps)
    end
end

# 一个函数 f 在 x 点的导数
"一个函数f在x点的导数"
function derivative(f, x, eps=0.02)
    return (f(x + eps) - f(x - eps)) / 2 / eps
end

# https://zenn.dev/cl17/articles/aec89318f89505

struct Point2
    x
    y
end

struct Interval
    left
    right
end
function x_diameter(i::Interval)
    i.right.x - i.left.x
end

"""
shubert_piyavskii 方法
## 参数
- f(x) 一元函数
- l 区间左值
- r 区间右值
- eps 误差限
## 返回值
Solution 结构体,含有
- point 最小值点
- value 最小值
"""
function shubert_piyavskii(f, l, r, eps)
    if l > r
        l, r = r, l
    end
    if r - l < eps
        p = (l + r) / 2
        return Solution(p, f(p))
    end
    # 凸函数的导数是递增的,导数最大值一定在两个端点
    max_derivative = max(abs(derivative(f, l)), abs(derivative(f, r)))
    k = max_derivative
    function possible_minimum(p1, p2)
        d = p2.x - p1.x
        # 有 dk - 2kx = y_2 - y_1
        min_point_delta = d / 2 + (p1.y - p2.y) / 2 / k
        min_point_y = p1.y - k * min_point_delta
        return min_point_y
    end
    mid = (l + r) / 2
    intervals = begin
        f_mid = f(mid)
        p_l = Point2(l, f(l))
        p_mid = Point2(mid, f_mid)
        p_r = Point2(r, f(r))
        i1 = Interval(p_l, p_mid)
        i2 = Interval(p_mid, p_r)
        PriorityQueue(i1 => possible_minimum(p_l, p_mid), i2 => possible_minimum(p_mid, p_r))
    end
    while true
        interval = dequeue!(intervals)
        mid_x = (interval.left.x + interval.right.x) / 2
        if x_diameter(interval) < eps
            return Solution(mid_x, f(mid_x))
        end
        mid_y = f(mid_x)
        p_l, p_r = interval.left, interval.right
        p_mid = Point2(mid_x, mid_y)
        i1 = Interval(p_l, p_mid)
        i2 = Interval(p_mid, p_r)
        push!(intervals, i1 => possible_minimum(p_l, p_mid))
        push!(intervals, i2 => possible_minimum(p_mid, p_r))
    end
end

end