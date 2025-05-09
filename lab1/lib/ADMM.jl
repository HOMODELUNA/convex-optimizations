module ADMM
export admm, Problem

using LinearAlgebra
include("NoGradientMethods.jl")
import .NoGradientMethods
"""
- x Px1
- z Qx1
## 项
- f(x) -> R 
- g(z) -> R
- A : MxP
- B : MxQ
- c : m 维向量
"""
struct Problem
    f
    g
    A
    B
    c
end
# https://zhuanlan.zhihu.com/p/448289351
"""
Admm 方法
## 参数
- problem 问题参数
- x0 初始x
- z0 初始z
"""
function admm(problem::Problem, x0, z0; rho=1.0, tol=1e-6, max_iter=1000)
    (; f, g, A, B, c) = problem
    x = x0
    z = z0
    (M, Q) = size(B)

    u = zeros(M)  # Dual variable

    for k = 1:max_iter
        # Step 1: Update x
        x = argmin_x(f, g, A, B, c, z, u, rho)

        # Step 2: Update z
        z_old = copy(z)
        z = argmin_z(g, B, A * x + u, c, rho)

        # Step 3: Update dual variable u
        u += A * x + B * z - c

        # Check convergence
        if norm(z - z_old) < tol && norm(A * x + B * z - c) < tol
            println("Converged after $k iterations.")
            break
        end
    end

    return x, z
end

function argmin_x(f, g, A, B, c, z, u, rho)
    L(x) = f(x) + g(A * x + B * z + u) + (rho / 2) * norm(A * x + B * z + u - c)^2
    result = optimize(L, zeros(size(A, 2)))
    return result
end

function argmin_z(g, B, Ax_plus_Bz_plus_u, c, rho)
    L(z) = g(z) + (rho / 2) * norm(Ax_plus_Bz_plus_u - c)^2
    result = optimize(L, zeros(size(B, 2)))
    return result
end

optimize = NoGradientMethods.hooke_jeeves

end