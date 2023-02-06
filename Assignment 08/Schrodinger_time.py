'''
Solve the time-dependent Schrodinger equation using both the Crank–Nicolson scheme and stable explicit scheme.
Consider the one-dimensional case and test it by applying it to the problem of a square well with a Gaussian initial state coming in from the left.
'''
import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

m = 1
hbar = 1
x_start = -20  # 设置研究波函数的空间起点与终点
x_end = 20
t_start = 0  # 设置研究波函数的时间起点与终点
t_end = 10
h_t = 0.01 #设置时间迭代的步长
h_x = 0.2  #设置空间迭代的步长
alpha = h_t/h_x**2

# 定义方势阱
def V(x):
    if -2 < x < 2:
        return -0.2
    else:
        return 0

# 定义高斯函数，其中k0=1，kesi0=-6，即高斯函数中点位于x=-6处
def Gauss(x):
    return math.sqrt(1/math.pi)*cmath.exp(1j*x-(x+6)**2/2)

T = np.arange(t_start, t_end+h_t, h_t)  # 生成t的坐标点
X = np.arange(x_start, x_end+h_x, h_x)  # 生成x的坐标点

######## Crank–Nicolson Scheme ########
# 生成初始时刻的波函数,以及波函数对应的概率密度分布
Psi1 = np.zeros((len(X), len(T)))*(1+0j)  #生成存储波函数的矩阵
Proba1 = np.zeros((len(X), len(T)))  #生成存储概率密度的矩阵
# 设置波函数的初态
psi = np.zeros((len(X), 1), dtype=complex)
for i in range(len(X)):
    psi[i][0] = Gauss(X[i])
# 定义矩阵iI+B:
Mat1 = np.zeros((len(X), len(X)), dtype=complex)
for i in range(len(X)):
    for j in range(len(X)):
        if i == j:
            Mat1[i, j] = 1j-V(X[i])*h_t/(2*m)-alpha*hbar/(2*m)
        if i == j+1 or i == j-1:
            Mat1[i, j] = alpha*hbar/(4*m)
# 定义矩阵iI-B:
Mat2 = np.zeros((len(X), len(X)), dtype=complex)
for i in range(len(X)):
    for j in range(len(X)):
        if i == j:
            Mat2[i, j] = 1j+V(X[i])*h_t/(2*hbar)+alpha*hbar/(2*m)
        if i == j+1 or i == j-1:
            Mat2[i, j] = -alpha*hbar/(4*m)
# 进行递推计算
Psi1[0:len(X), 0:1] = psi
Proba1[0:len(X), 0:1] = (abs(psi))**2
for i in range(1, len(T)):
    tem = np.dot(Mat2, psi)
    psi = np.dot(np.linalg.inv(Mat1), tem)
    Psi1[0:len(X), i:i+1] = psi  #求出下一时刻波函数的值
    Proba1[0:len(X), i:i+1] = (abs(psi))**2  #求出下一时刻概率密度的值
# 绘制概率密度随时间的演化图
figure = plt.figure()
ax = figure.add_subplot(111, projection='3d')
t, x = np.meshgrid(T, X)
ax.plot_surface(t, x, Proba1, rstride=1, cstride=1, cmap='rainbow')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('Probability')
plt.show()


######## Stable Explicit Scheme ########
# 生成初始时刻的波函数,以及波函数对应的概率密度分布
Psi2 = np.zeros((len(X), len(T)))*(1+0j)  #生成存储波函数的矩阵
Proba2 = np.zeros((len(X), len(T)))  #生成存储概率密度的矩阵
# 设置波函数最初的两个态，初态已知，第一个态由Crank-Nicolson Scheme计算
psi = np.zeros((len(X), 1), dtype=complex)
psi_1 = np.zeros((len(X), 1), dtype=complex)
for i in range(len(X)):
    psi_1[i][0]=Psi1[i][0]
    psi[i][0]=Psi1[i][1]
    Psi2[i][0] = psi_1[i][0]
    Proba2[i][0] = (abs(psi_1[i][0]))**2
    Psi2[i][1] = psi[i][0]
    Proba2[i][1] = (abs(psi[i][0]))**2
# 设置矩阵A
Mat3 = np.zeros((len(X), len(X)), dtype=complex)
for i in range(len(X)):
    for j in range(len(X)):
        if i == j:
            Mat3[i, j] = -2j*alpha*hbar**2/m-2j*V(X[i])*h_t
        if i == j+1 or i == j-1:
            Mat3[i, j] = 1j*alpha*hbar**2/m
# 开始递推计算
for i in range(2,len(T)): 
    tem=psi
    psi = psi_1+np.dot(Mat3,psi)
    psi_1=tem
    Psi2[0:len(X), i:i+1] = psi
    Proba2[0:len(X), i:i+1] = (abs(psi))**2
# 绘制概率密度随时间的演化图
figure = plt.figure()
ax = figure.add_subplot(111, projection='3d')
t, x = np.meshgrid(T, X)
ax.plot_surface(t, x, Proba2, rstride=1, cstride=1, cmap='rainbow')
ax.set_xlabel('t')
ax.set_ylabel('x')
ax.set_zlabel('Probability')
plt.show()
