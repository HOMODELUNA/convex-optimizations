
function f(x)
    x[1]^2 + x[2]^2 - x[1] * x[2] - 10 * x[1] - 4 * x[2] + 60
end

function grad_f(x)
    [2 * x[1] - x[2] - 10, 2 * x[2] - x[1] - 4]
end

include("lib/davidon_flether_powell.jl")
# BFGS optimization algorithm
function bfgs(f, gradient, x0; eps=1e-6, max_iter=1000)
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

struct CGArgs
    A
    b
end

function apply(args::CGArgs, x)
    0.5 * transpose(x) * args.A * x - dot(args.b, x)
end
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
        return iter(k - 1, g_next, p_next, x_next)
    end
    return iter(iter_max, grad, p0, x0)

end
x_0 = [0, 0]

@show davidon_flether_powell(f, grad_f, x_0)
@show bfgs(f, grad_f, x_0)


args = begin
    A = [2 -1; -1 2]
    b = [10, 4]
    CGArgs(A, b)
end


@show conjugate_gradient(args, x_0)

