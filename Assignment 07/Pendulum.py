'''
Write a code to numerically solves the motion of a simple pendulum using Euler’s method, 
midpoint method, RK4, Euler-trapezoidal method (implement these methods by yourself). 
Plot the angle and total energy as a function of time. Explain the results.
'''

import numpy as np
import matplotlib.pyplot as plt
import math

h=0.01  #步长
l=1  #绳长
m=1  #小球质量
g=9.8  #重力加速度
print('Please input the times you want to calculate:')
N=int(input())  #总共计算的次数

x=[i*h for i in range(N+1)]  #时间

def Energy(y1,y2):
    E=[m*l**2/2*y2[i]**2+m*g*l*(1-math.cos(y1[i])) for i in range(N+1)]
    return E

y1_Euler=np.zeros(N+1)  #角度
y2_Euler=np.zeros(N+1)  #角速度
y1_Euler[0]=0.5  #初始角度设为0.5rad
y2_Euler[0]=0  #小球由静止释放
E_Euler=np.zeros(N+1)
def Euler(N,y1,y2):  #Euler's Method
    for i in range(N):
        y1[i+1]=y1[i]+h*y2[i]
        y2[i+1]=y2[i]-(h*g/l)*math.sin(y1[i])

y1_Mid=np.zeros(N+1)  #角度
y2_Mid=np.zeros(N+1)  #角速度
y1_Mid[0]=0.5  #初始角度设为0.5rad
y2_Mid[0]=0  #小球由静止释放
E_Mid=np.zeros(N+1)
def Midpoint(N,y1,y2):  #Midpoint Method
    for i in range(N):
        delta_y1=h*y2[i]
        delta_y2=-h*g/l*math.sin(y1[i])
        y1[i+1]=y1[i]+h*(y2[i]+delta_y2/2)
        y2[i+1]=y2[i]-h*g/l*math.sin(y1[i]+delta_y1/2)

y1_Rk4=np.zeros(N+1)  #角度
y2_Rk4=np.zeros(N+1)  #角速度
y1_Rk4[0]=0.5  #初始角度设为0.5rad
y2_Rk4[0]=0  #小球由静止释放
E_Rk4=np.zeros(N+1)
def Rk4(N,y1,y2):  #4th order Runge-Kutta Method
    for i in range(N):
        S1_y1=y2[i]  #求出初始位置的斜率，S1
        S1_y2=-g/l*math.sin(y1[i])
        y1_mid=y1[i]+h/2*S1_y1  #由S1求出中点函数值
        y2_mid=y2[i]+h/2*S1_y2
        S2_y1=y2_mid  #求出中点的斜率，S2
        S2_y2=-g/l*math.sin(y1_mid)
        y1_mid=y1[i]+h/2*S2_y1  #由中点斜率（S2）重新计算中点函数值
        y2_mid=y2[i]+h/2*S2_y2
        S3_y1=y2_mid
        S3_y2=-g/l*math.sin(y1_mid)  #再重新计算中点斜率，S3
        y1[i+1]=y1[i]+h*S3_y1  #计算末端函数值
        y2[i+1]=y2[i]+h*S3_y2
        S4_y1=y2[i+1]  #计算末端斜率，S4
        S4_y2=-g/l*math.sin(y1[i+1])
        y1[i+1]=y1[i]+h/6*(S1_y1+2*S2_y1+2*S3_y1+S4_y1)
        y2[i+1]=y2[i]+h/6*(S1_y2+2*S2_y2+2*S3_y2+S4_y2)

y1_EulerTra=np.zeros(N+1)  #角度
y2_EulerTra=np.zeros(N+1)  #角速度
y1_EulerTra[0]=0.5  #初始角度设为0.5rad
y2_EulerTra[0]=0  #小球由静止释放
E_EulerTra=np.zeros(N+1)
def Euler_Trape(N,y1,y2):  #Euler-trapezoidal Method
    for i in range(N):
        y1_1=y1[i]  #y1_1与y2_2表示本次估算值
        y2_1=y2[i]
        y1_2=y1[i]+h*y2[i]  #y1_2与y2_2表示下一次估算值
        y2_2=y2[i]-h*g/l*math.sin(y1[i])
        while abs(y1_1-y1_2)>1e-5 or abs(y2_1-y2_2)>1e-5:  #如果两次估算值差距较大，则继续迭代
            y1_1=y1_2
            y2_1=y2_2
            y1_2=y1[i]+h/2*(y2[i]+y2_1)
            y2_2=y2[i]-h*g/2/l*(math.sin(y1[i])+math.sin(y1_1))
        y1[i+1]=y1_2
        y2[i+1]=y2_2


Euler(N,y1_Euler,y2_Euler)
Midpoint(N,y1_Mid,y2_Mid)
Rk4(N,y1_Rk4,y2_Rk4)
Euler_Trape(N,y1_EulerTra,y2_EulerTra)

E_Euler=Energy(y1_Euler,y2_Euler)
E_Mid=Energy(y1_Mid,y2_Mid)
E_Rk4=Energy(y1_Rk4,y2_Rk4)
E_EulerTra=Energy(y1_EulerTra,y2_EulerTra)

fig1,ax1=plt.subplots()
ax1.set_xlabel('t')
ax1.set_ylabel('angle')
ax1.set_title('angle-t graph of the pendulum')
ax1.plot(x,y1_Euler, color='b', label="Euler")
ax1.plot(x,y1_Mid, color='r', label="Midpoint")
ax1.plot(x,y1_Rk4, color='y', label="Rk4")
ax1.plot(x,y1_EulerTra, color='g', label="Euler Trapezoid")
ax1.legend(loc='best')
plt.show()

fig2,ax2=plt.subplots()
ax2.set_xlabel('t')
ax2.set_ylabel('energy')
ax2.set_title('energy-t graph of the pendulum')
ax2.plot(x,E_Euler, color='b', label="Euler")
ax2.legend(loc='best')
plt.show()

fig3,ax3=plt.subplots()
ax3.set_xlabel('t')
ax3.set_ylabel('energy')
ax3.set_title('energy-t graph of the pendulum')
ax3.plot(x,E_Mid, color='r', label="Midpoint")
ax3.legend(loc='best')
plt.show()

fig4,ax4=plt.subplots()
ax4.set_xlabel('t')
ax4.set_ylabel('energy')
ax4.set_title('energy-t graph of the pendulum')
ax4.plot(x,E_Rk4, color='y', label="Rk4")
ax4.legend(loc='best')
plt.show()

fig5,ax5=plt.subplots()
ax5.set_xlabel('t')
ax5.set_ylabel('energy')
ax5.set_title('energy-t graph of the pendulum')
ax5.plot(x,E_EulerTra, color='g', label="Euler Trapezoid")
ax5.legend(loc='best')
plt.show()