"""
'encoding'表示每个变量的编码方式，它的值可以是一个标量或行向量，
且每维的值可以为1（实数）、2（整数）、3（标签）、4（二进制数）或5（序
列编号）。算法针对不同的编码方式可能使用不同的算子来产生解。
"""
@enum Encoding Real Integer Tag Binary SerialNumber

"""
## 目标函数
``f(x) = (f_1(x),f_2(x),...,f_m(x))``
## 结构体的项
- targets 一个目标函数数组,每项为``f(x)``的格式,其中x为``(x_1,x_2,...,x_m)``
"""
struct Problem
    targets
    constraints
end

p = Problem(1,1)
p.constraints