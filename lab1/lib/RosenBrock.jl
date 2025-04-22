module RosenBrock
a = 1
b = 100
function f(xs)
    x = xs[1]
    y = xs[2]
    return (a - x)^2 + b * (y - x^2)^2
end

function g(x)
    x_grad = -2 * (a - x[1]) - 4 * b * x[1] * (x[2] - x[1]^2)
    y_grad = 2 * b * (x[2] - x[1]^2)
    return [x_grad, y_grad]
end

"""
RosenBrock 函数的二阶导数
## 输入
- x 向量
## 返回值
- 一个 2x2 矩阵
"""
function g2(x)
    xx = 2 - 4*b*x[2] + 12 * b * x[1]^2
    xy = -4 * b * x[1]
    yx = - 4 * b * x[1]
    yy = 2 * b
    return [xx xy; yx yy]
end

export f, g, g2

# 最小值点 [a, a^2]
end