include("types.jl")

# https://blog.csdn.net/tianhai12/article/details/130491491
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