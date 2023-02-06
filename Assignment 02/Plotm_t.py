'''
m(t)=tanh(m(t)/t)
for a given t, solve m, plot m as a function of t
'''
'''
  1.首先研究方程f(m,t)=tanh(m/t)-m解的情况,即以0.5为步长取定t分别为[-2,2]中的数,作出f(m)的图像
  2.由图可以知道,当t≤0与t≥1的时候,函数都只有m=0这一个零点
    当0<t<1时,函数有三个零点,除m=0这一平凡解之外,还有关于原点对称的两个零点,并且这两个零点的绝对值总在0到1之间
  3.考虑真实的物理场景:约化温度t>0,m(t)的两个非零解大小相等，正负代表磁化方向,m(t)的平凡解则不予考虑
  4.因此我们只需要着重分析当t∈(0,1)时m∈(0,1)的解,这一部分的解可以由二分法求得
'''
import numpy as np
import matplotlib.pyplot as plt
import math


def f(m, t):  # 定义函数
    return (math.tanh(m/t)-m)


def Plot_f(t):   #对于一个确定的t，绘制出f(m,t)的函数图像，以探究m与t是否一一对应
    x=np.arange(-5,5,0.005)  #横坐标的值从-5取到5，间隔为0.005
    y=[]
    for m in x:
        y.append(f(m,t))   #将对应的f(m)的值放入y中
    plt.xlabel('m')  #横坐标为m
    plt.ylabel('f(m)')   #纵坐标为f(m)
    plt.plot(x,y)  #画出图像
    plt.ylim(-2,2)
    plt.show()


def FindRoot(a, b, t):  #对于给定的t的值，利用二分法求解对应的m的值
    if t>1:   #由之前的分析可知，当t>1时，只有平凡解m=0
        return 0
    else:
        while abs(a-b) > 0.0001:
            if f(a, t) * f((a+b)/2, t) < 0:
                b = (a+b)/2
            else:
                a = (a+b)/2
        return a


t_list=np.linspace(0.000001,1.5,500)  #取因变量t的取值，从0到1.5，一共取为500个点。为防止出现分母为0，将t的初始值取为0.000001
m_list=[]
for t in t_list:
    m_list.append(FindRoot(0.000001,1,t))   #由t的取值用二分法求解出对应的m的值
                                            #注意此处区间从0.000001的原因是，当a=0时，将不能取到新的区间
plt.plot(t_list,m_list)
plt.ylim(0,1)
plt.xlabel('t')
plt.ylabel('m(t)')
plt.show()