using LinearAlgebra

f(x) = x[1] - x[2] + 2 * x[1]^2 + 2 * x[1] * x[2] + x[2]^2

struct CGArgs
    A
    b
end

function apply(args::CGArgs, x)
    0.5 * transpose(x) * args.A * x - dot(args.b, x)
end

"""
## 注意 args.A 应当为对称的
"""
function cg_gradient(args::CGArgs, x)
    (; A, b) = args
    return A * x - b
end

Arg = begin
    A = [4 2; 2 2]
    b = [-1, 1]
    CGArgs(A, b)
end
X_INIT = [0, 0]
# https://zhuanlan.zhihu.com/p/178461470 

# https://zhuanlan.zhihu.com/p/234950550

#  这一参考最好 https://zhuanlan.zhihu.com/p/23811968

struct Solution
    point
    value
end

"""
``ϕ(x) = \\frac{1}{2} x^TAx - x^Tb``
"""
function conjugate_gradient(args::CGArgs, x0; eps=1e-6, iter_max=100)
    (; A, b) = args
    grad = A * x0 - b
    p0 = -grad
    function iter(k, g_k, pk, x)
        if k <= 0
            return Solution(x, apply(args, x))
        end
        if norm(pk, 1) < eps
            return Solution(x, apply(args, x))
        end
        grad_square = dot(g_k, g_k)
        alpha_k = grad_square / (transpose(pk) * A * pk)
        x_next = x + alpha_k * pk
        if norm(x_next - x) < eps
            return Solution(x_next, apply(args, x_next))
        end
        g_next = g_k + alpha_k * A * pk
        beta_k = dot(g_next, g_next) / grad_square
        p_next = -g_next + beta_k * pk
        return iter(k-1,g_next,p_next,x_next)
    end
    return iter(iter_max, grad, p0, x0)

end

@show apply(Arg, X_INIT)
@show conjugate_gradient(Arg, X_INIT)