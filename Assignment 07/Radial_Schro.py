'''
Write a code to numerically solves radial Schrödinger equation for:
1. V(r)=-1/r
2. V(r)=V_loc
'''

import numpy as np
import matplotlib.pyplot as plt
import math


N=1000  #设置求解区间的精度
r_max=35  #设置最大的半径值
r_min=1e-10
t_max=6  #设置求解区间的长度
r_p=(r_max-r_min)/(math.exp(t_max)-1)
h=t_max/N
t=np.linspace(0,t_max,N+1)

def r(t):
    return r_p*(math.exp(t)-1)+r_min

def V_hydrogen(t):
        return -1/r(t)

def V_Li(t):
    Z_ion=3
    r_loc=0.4
    C1=-14.0093922
    C2=9.5099073
    C3=-1.7532723
    C4=0.0834586
    return -Z_ion/r(t)*math.erf(r(t)/math.sqrt(2)/r_loc)+math.exp(-0.5*(r(t)/r_loc)**2)*(C1+C2*(r(t)/r_loc)**2+C3*(r(t)/r_loc)**4+C4*(r(t)/r_loc)**6)

def RK4(V,E,y1,y2,l):
    def deri_y2(y1,y2,t):
        return y2-(math.exp(2*t))*2*r_p**2*(E-V(t)-l*(l+1)/(2*r(t)**2))*y1
    for i in range(N,0,-1):
        S1_y1=y2[i]
        S1_y2=deri_y2(y1[i],y2[i],t[i])
        y1_mid=y1[i]-h/2*S1_y1  #由S1求出中点函数值
        y2_mid=y2[i]-h/2*S1_y2
        S2_y1=y2_mid  #求出中点的斜率，S2
        S2_y2=deri_y2(y1_mid,y2_mid,t[i])
        y1_mid=y1[i]-h/2*S2_y1  #由中点斜率（S2）重新计算中点函数值
        y2_mid=y2[i]-h/2*S2_y2
        S3_y1=y2_mid
        S3_y2=deri_y2(y1_mid,y2_mid,t[i])  #再重新计算中点斜率，S3
        y1[i-1]=y1[i]-h*S3_y1  #计算末端函数值
        y2[i-1]=y2[i]-h*S3_y2
        S4_y1=y2[i-1]  #计算末端斜率，S4
        S4_y2=deri_y2(y1[i-1],y2[i-1],t[i-1])
        y1[i-1]=y1[i]-h/6*(S1_y1+2*S2_y1+2*S3_y1+S4_y1)
        y2[i-1]=y2[i]-h/6*(S1_y2+2*S2_y2+2*S3_y2+S4_y2)
    return y1[0]

def shoot(V,E1,l):
    y1=np.zeros(N+1)
    y2=np.zeros(N+1)
    y1[N]=0  #r_max处的波函数设为0
    y2[N]=1.0  #初始位置的导数值并不影响最终结果
    E2=E1-1e-3
    y1_first=RK4(V,E1,y1,y2,l)
    while abs(y1_first)>1e-4:
        y1_first_=RK4(V,E2,y1,y2,l)
        tmp = E2
        E2=E2-y1_first_*(E2-E1)/(y1_first_-y1_first)
        E1=tmp
        y1_first=y1_first_
    return E2,y1
'''
#绘制E~u0(E)的关系图，从而初步判断能量范围
E_list=np.linspace(-0.5,0,100)
y1_list=[]
y1_test=np.zeros(N+1)
y2_test=np.zeros(N+1)
y1_test[N]=0  #r_max处的波函数设为0
y2_test[N]=1.0  #初始位置的导数值并不影响最终结果
for i in range(len(E_list)):
    y1_list.append(RK4(V_Li,E_list[i],y1_test,y2_test,1))  #RK4的第一个参数代表势能函数，可选V_hydrogen或V_Li
plt.plot(E_list,y1_list)                                         #最后一个参数代表了l，可选l=0或l=1
plt.ylim(-1e2,1e2)
plt.show()
'''
###求解氢原子波函数
#求解n=1,l=0的氢原子1s态
E_hydorgen1,y1_hydrogen1=shoot(V_hydrogen,-0.8,0)
#求解n=2,l=0的氢原子2s态
E_hydorgen2,y1_hydrogen2=shoot(V_hydrogen,-0.3,0)
#求解n=2,l=1的氢原子2p态
E_hydorgen3,y1_hydrogen3=shoot(V_hydrogen,-0.3,1)
#输出最小的三个能量值
print(E_hydorgen1,E_hydorgen2,E_hydorgen3)
#绘制波函数图像，红色为R(r)图像，蓝色为u(r)图像
r_list=[r(t[i]) for i in range(N+1)]
R_hydrogen1=[-y1_hydrogen1[i]/r_list[i] for i in range(N+1)]
R_hydrogen2=[y1_hydrogen2[i]/r_list[i] for i in range(N+1)]
R_hydrogen3=[-y1_hydrogen3[i]/r_list[i] for i in range(N+1)]
r_list.pop(0)
R_hydrogen1.pop(0)
R_hydrogen2.pop(0)
R_hydrogen3.pop(0)
#绘制第一幅图
fig1,ax1=plt.subplots()
ax1.set_xlabel('r')
ax1.set_ylabel('R(r)')
ax1.plot(r_list,R_hydrogen1,color='r',label="R(r)")
ax1.plot(r_list,-y1_hydrogen1[1:],color='b',label="u(r)")
ax1.set_title('Hydrogen n=1, l=0, 1s')
ax1.legend(loc='best')
plt.show()
#绘制第二幅图
fig2,ax2=plt.subplots()
ax2.set_xlabel('r')
ax2.set_ylabel('R(r)')
ax2.plot(r_list,R_hydrogen2,color='r',label="R(r)")
ax2.plot(r_list,y1_hydrogen2[1:],color='b',label="u(r)")
ax2.set_title('Hydrogen n=2, l=0, 2s')
ax2.legend(loc='best')
plt.show()
#绘制第三幅图
fig3,ax3=plt.subplots()
ax3.set_xlabel('r')
ax3.set_ylabel('R(r)')
ax3.plot(r_list,R_hydrogen3,color='r',label="R(r)")
ax3.plot(r_list,-y1_hydrogen3[1:],color='b',label="R(r)")
ax3.set_title('Hydrogen n=2, l=1, 2p')
ax3.legend(loc='best')
plt.show()

###求解锂原子波函数
#求解l=0的锂原子态
E_Li1,y1_Li1=shoot(V_Li,-1,0)

#求解l=0的锂原子态
E_Li2,y1_Li2=shoot(V_Li,-0.4,0)
#求解l=1的锂原子态
E_Li3,y1_Li3=shoot(V_Li,-1,1)
#输出最小的三个能量值
print(E_Li1,E_Li2,E_Li3)
#绘制波函数图像，红色为R(r)图像，蓝色为u(r)图像
r_list=[r(t[i]) for i in range(N+1)]
R_Li1=[-y1_Li1[i]/r_list[i] for i in range(N+1)]
R_Li2=[y1_Li2[i]/r_list[i] for i in range(N+1)]
R_Li3=[y1_Li3[i]/r_list[i] for i in range(N+1)]
r_list.pop(0)
R_Li1.pop(0)
R_Li2.pop(0)
R_Li3.pop(0)

#绘制第一幅图
fig1,ax1=plt.subplots()
ax1.set_xlabel('r')
ax1.set_ylabel('R(r)')
ax1.plot(r_list,R_Li1,color='r',label="R(r)")
ax1.plot(r_list,-y1_Li1[1:],color='b',label="u(r)")
ax1.set_title('Li, l=0, E_Li1')
ax1.legend(loc='best')
plt.show()
#绘制第二幅图
fig2,ax2=plt.subplots()
ax2.set_xlabel('r')
ax2.set_ylabel('R(r)')
ax2.plot(r_list,R_Li2,color='r',label="R(r)")
ax2.plot(r_list,y1_Li2[1:],color='b',label="u(r)")
ax2.set_title('Li, l=0, E_Li2')
ax2.legend(loc='best')
plt.show()
#绘制第三幅图
fig3,ax3=plt.subplots()
ax3.set_xlabel('r')
ax3.set_ylabel('R(r)')
ax3.plot(r_list,R_Li3,color='r',label="R(r)")
ax3.plot(r_list,y1_Li3[1:],color='b',label="R(r)")
ax3.set_title('Li, l=1, E_Li3')
ax3.legend(loc='best')
plt.show()

