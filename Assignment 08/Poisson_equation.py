'''
Consider the Poisson equation from electrostatics on a rectangular geometry with x ∈ [0, Lx ] and y ∈[0, Ly ]. 
Write a program that solves this equation using the relaxation method. Test your program with: 
(a) ρ(x, y) = 0, φ(0, y) = φ(Lx , y) =φ(x, 0) = 0, φ(x, Ly ) = 1V, Lx = 1m, and Ly = 1.5m
(b) ρ(x, y)/ɛ0 = 1 V/m^2, φ(0, y) = φ(Lx , y) = φ(x, 0) = φ(x, Ly ) = 0, and Lx = Ly = 1 m
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

h=0.02  #设置求解的步长

#定义迭代函数
def SOR(V,w,X,Y,h,p):
    for k in range(3000):  #如果迭代超过最大迭代次数，则停止迭代
        signal=0  #用于表示迭代前后的差距
        for i in range(1,len(Y)-1):  #遍历每一行
            for j in range(1,len(X)-1):  #遍历每一列
                tem=V[i][j]
                V[i][j]=(1-w)*V[i][j]+w/4*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]+h**2*p)
                signal+=abs(tem-V[i][j])
        signal=signal**0.5
        if signal<1e-5:  #如果迭代前后的差距在误差范围内，则推出迭代，返回结果
                break
    return V

#定义初始化函数，用于设定求解问题的边界条件和电荷分布
def Poisson(Lx,Ly,bx1,bx2,by1,by2,p):   #参数Lx，Ly为区间长度；bx1，bx2为x的两个边界值；by1，by2为y的两个边界值；p为电荷密度
    w=2.0/(1+np.pi/Lx)  #定义SOR算法中的over-relaxation parameter
    X=np.arange(0,Lx+h,h)  #生成x方向的坐标点
    Y=np.arange(0,Ly+h,h)  #生成y方向的坐标点
    V=np.zeros((len(Y),len(X)))  #生成电势的格点
    # 设定初始条件
    for i in range(len(Y)):  #x=0与x=Lx处的初始条件
        V[i][len(X)-1] = bx2
        V[i][0] = bx1
    for i in range(len(X)):  #y=0与y=Ly处的初始条件
        V[0][i] = by1
        V[len(Y)-1][i] = by2
    for i in range(1,len(Y)-1):  #在边界以内的各点随机生成一个试解
        for j in range(1,len(X)-1):
            V[i][j]=np.random.uniform(0,1)
    V=SOR(V,w,X,Y,h,p)  #将上述各值输入迭代函数中进行迭代求解
    #绘制势能3D图
    figure = plt.figure()
    ax = figure.add_subplot(111, projection='3d')
    x,y=np.meshgrid(X,Y)
    ax.plot_surface(x,y,V,rstride=1,cstride=1,cmap='rainbow')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('V')
    plt.show()
    #绘制势能平面投影图
    CS = plt.contourf(X,Y,V, linewidth=2, cmap=mpl.cm.jet)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar(CS)
    plt.show()

#求解第一个问题
Poisson(1,1.5,0,0,0,1,0)
#求解第二个问题
Poisson(1,1,0,0,0,0,1)