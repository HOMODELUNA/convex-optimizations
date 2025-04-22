# 实验1

## 实验要求


## 测试函数

- [Booth 函数](lib/Booth.jl)
- [RosenBrock 函数](lib/RosenBrock.jl)

## 无导数求解算法

见 [lib/NoGradientMethods.jl](lib/NoGradientMethods.jl)
实现方法:
- 轮换坐标法
- Hooke Jeeves 搜索法

测试方法:
```sh
julia use-no-gradient.jl
```
## 梯度算法
见 [lib/GradientMethods.jl](lib/GradientMethods.jl)

实现方法:
- 梯度下降法
- 梯度下降法(带动量)
- 共轭梯度法
- 共轭梯度法(带动量)


测试方法:
```sh
julia use-gradient-descend.jl
```
## 牛顿法

见 [lib/NewtonMethods.jl](lib/NewtonMethods.jl)

实现了
- 牛顿法
- 阻尼牛顿法
- [] 三次正则化牛顿法

测试方法:
```sh
julia use-newton.jl
```