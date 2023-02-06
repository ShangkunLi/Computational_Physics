'''
 Newton interpolation of 
 (i) 10 equal spacing points of cos(x) within [0,p], 
 (ii) 10 equal spacing points 1/(1+25x2) within [-1,1]. 
 Compare the results with the cubic spline interpolation.
'''
'''
  本题由主要由两个部分完成整体功能:
  Newton Interpolation完成牛顿插值法的功能
  Cubic Spline interpolation完成三阶样条插值法的功能
  1. Newton Interpolation
  (1) 首先定义Bracket_fun()函数，用于求出牛顿插值法中的Bracket Function，即f[xi1,xi2,...,xij]，本题中采用递归算法求得
  (2) 定义Newton_Intert()函数，可以返回某一自变量值下对应插值函数的值
  (3) 定义Plot_Newton()函数，用于绘制计算得到的函数图像
  2. Cubic Spline interpolation
  (1) 首先定义Second_Derivative()函数，通过求解一个三对角矩阵，得到样本点处的二阶导函数值
  (2) 定义Plot_Spline()函数，绘制出计算得到的函数图像
'''


import numpy as np
import matplotlib.pyplot as plt
import math
N = 10  #设置用于差值的样本点个数
#生成函数cos(x)的差值点
x1_input = np.linspace(0, math.pi, N)
f1_input = [math.cos(i) for i in np.linspace(0, math.pi, N)]
#生成函数1/(1+25*x^2)的差值点
x2_input = np.linspace(-1, 1, N)
f2_input = [1.0/(1.0+25.0*i**2) for i in np.linspace(-1, 1, N)]


# ------------------Newton interpolation------------------

def Bracket_fun(xi_input, fi_input):  #定义Bracket_fun()函数，利用递归求出Bracket Function
    if len(xi_input) > 2 and len(fi_input) > 2:  #如果数组中元素数量大于2，则使用Bracket Function之间的递推关系进行递归降阶
        return (Bracket_fun(xi_input[1:len(xi_input)], fi_input[1:len(fi_input)])-Bracket_fun(xi_input[:len(xi_input)-1], fi_input[:len(xi_input)-1]))/float(xi_input[-1]-xi_input[0])
    else:  #如果数组中仅剩两个元素，则通过Bracket Function的定义计算
        return (fi_input[1]-fi_input[0])/float(xi_input[1]-xi_input[0])


def Newton_Inter(xi_input, fi_input, x):  #定义Newton_Inter()，计算对应x点处的函数值
    r = 1.0
    f = fi_input[0]
    for i in range(N):
        r *= (x-xi_input[i])  #计算牛顿差值公式中的多项式部分
        f += Bracket_fun(xi_input[0:i+2], fi_input[0:i+2])*r  #在各个多项式前面乘以对应的系数，即Bracket Function
    return f


def Plot_Newton(xi_input, fi_input, xi_newton):  #定义  Plot_Newton()，绘制函数图像
    fi_newton = []
    for x in xi_newton:
        fi_newton.append(Newton_Inter(xi_input, fi_input, x))
    plt.plot(xi_newton, fi_newton, color='r')

# --------------------------------------------------------


# ----------------Cubic Spline interpolation---------------

def Second_Derivative(xi_input, fi_input):  #定义Second_Derivative()函数，求解样本点处的二阶导函数值
    Solution = np.zeros((N, 1))  #声明一个10*1的列向量，用于存储二阶导函数值
    Matrix = np.zeros((N, N))   #声明一个10*10的矩阵，是用于求解二阶导函数值的系数矩阵
    Vector = np.zeros((N, 1))   #声明一个10*1的列向量，是线性方程组Ax=b中的b
    for i in range(1, N-1):  #通过函数在样本点处一阶导函数连续对系数矩阵A和向量b进行赋值
        Matrix[i][i-1] = xi_input[i]-xi_input[i-1]
        Matrix[i][i] = 2*(xi_input[i+1]-xi_input[i-1])
        Matrix[i][i+1] = xi_input[i+1]-xi_input[i]
        Vector[i][0] = 6.0*(fi_input[i+1]-fi_input[i])/(xi_input[i+1]-xi_input[i]) + \
            6.0*(fi_input[i-1]-fi_input[i])/(xi_input[i]-xi_input[i-1])
    Matrix[0][0] = 1.0
    Matrix[N-1][N-1] = 1.0
    Vector[0][0] = 0.0
    Vector[N-1][0] = 0.0
    Solution = np.matmul(np.linalg.inv(Matrix), Vector)   #求解二阶导函数值
    return Solution


def Plot_Spline(Solution, xi_input, fi_input):      #定义Plot_Spline()，绘制函数图像
    fi_spline = []
    for i in range(N-1):  #对与不同段的x取值，用不同的函数进行求解
        xi_spline = np.linspace(xi_input[i], xi_input[i+1], 50)
        fi_spline = [Solution[i][0]*(xi_input[i+1]-x)**3/6/(xi_input[i+1]-xi_input[i])+Solution[i+1][0]*(-xi_input[i]+x)**3/6/(xi_input[i+1]-xi_input[i])
                     + (fi_input[i]/(xi_input[i+1]-xi_input[i])-Solution[i]
                        [0]*(xi_input[i+1]-xi_input[i])/6)*(xi_input[i+1]-x)
                     + (fi_input[i+1]/(xi_input[i+1]-xi_input[i])-Solution[i+1][0]*(xi_input[i+1]-xi_input[i])/6)*(x-xi_input[i]) for x in xi_spline]
        plt.plot(xi_spline, fi_spline,color='b')
# ---------------------------------------------------------

#绘制cos(x)的函数图像，红线为牛顿插值法的结果，蓝线为样条插值法的结果
x1_newton = np.linspace(0, math.pi, 200)
Plot_Newton(x1_input, f1_input, x1_newton)
Plot_Spline(Second_Derivative(x1_input, f1_input), x1_input, f1_input)
plt.scatter(x1_input, f1_input,color='g')
plt.title('cos(x)')
plt.show()

#绘制1/(1+25*x^2)的函数图像，红线为牛顿插值法的结果，蓝线为样条插值法的结果
x2_newton = np.linspace(-1, 1, 200)
Plot_Newton(x2_input, f2_input, x2_newton)
Plot_Spline(Second_Derivative(x2_input, f2_input), x2_input, f2_input)
plt.scatter(x2_input, f2_input,color='g')
plt.title('1/(1+25*x^2)')
plt.show()
