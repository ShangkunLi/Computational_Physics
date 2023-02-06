'''
Compute ∫(R_3s)^2*r^2*dr for Si atom (Z=14) with Simpson's rule using two different radial grids:
(1) Equal spacing grids: r[i]=(i-1)h; i = 1, ..., N (try different N)
(2) A nonuniform integration grid, more finely spaced at small r than at large r: r[i] = r0 (e^t[i]-1); t[i]=(i-1)h; i = 1, ..., N (One typically choose r0 = 0.0005 a.u., try different N).
(3) Find out which one is more efficient, and discuss the reason.
'''
'''
本题主要利用Simpson法计算数值积分
通过对自变量r进行不同的变量代换，我们得到不同的待积函数
对于两种不同的待积函数，我们分别取其积分函数点个数为N1和N2
'''
import numpy as np
import math
import matplotlib.pyplot as plt

print('Plese input two numbers of points:(must be odd number, eg. 107 101)')  #区间数必须为偶数，因此函数点数必须为奇数
N1,N2=map(int,input().split())
r0=0.0005
r=np.linspace(0,40,N1)  #将r分N1个点
t=np.linspace(0,math.log(1+40/r0),N2)  #将t分为N2个点
intergra1=0.0   #函数f1的积分结果
intergra2=0.0  #函数f2的积分结果
def f1(r):  #定义以r为积分变量的函数
    return 1.0/(81*3)*(6-56*r+(28*r/3)**2)*math.pow(14,1.5)*math.exp(-14*r/3)*r**2
def f2(t):  #定义以t为积分变量的函数
    return 1.0/(81*3)*(6-56*r0*(math.exp(t)-1)+(28*r0/3*(math.exp(t)-1))**2)*math.pow(14,1.5)*math.exp(-14*r0/3*(math.exp(t)-1))*r0**2*(math.exp(t)-1)**2*r0*math.exp(t)
for i in range(0,N1-2,2):  #利用Simpson公式进行积分
    intergra1+=(r[i+1]-r[i])/3*(f1(r[i])+4*f1(r[i+1])+f1(r[i+2]))
for i in range(0,N2-2,2):  #利用Simpson公式进行积分
    intergra2+=(t[i+1]-t[i])/3*(f2(t[i])+4*f2(t[i+1])+f2(t[i+2]))
print('intergra1=',intergra1)
print('intergra2=',intergra2)

'''
#这部分代码用于绘制f1(r)与f2(t)在积分区间内的函数图像
x1_axis=np.linspace(0,40,500)
x2_axis=np.linspace(0,math.log(1+40/r0),500)
y1_axis=[f1(i) for i in x1_axis]
y2_axis=[f2(i) for i in x2_axis]
plt.plot(x1_axis,y1_axis,color='r',label='f1(r)')
plt.plot(x2_axis,y2_axis,color='b',label='f2(t)')
plt.legend()
plt.show()
'''
