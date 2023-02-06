'''
Write a MC code for a 3D Face-Centered Cubic lattice using the Heisenberg spin model (adopt periodic boundary condition).
Estimate the ferromagnetic Curie temperature. 
'''
import numpy as np
from math import cos, sin,acos
import math
import random
import matplotlib.pyplot as plt

N=10 #设置格点的个数
T=10  #设置求解的最大温度
max_loop_times=N**3*50  #设置每一温度最大的迭代次数，以此求得每一温度下稳定的结果
axi_tem=np.linspace(0.01,T,40)  #设置求解的温度范围
energy=[]

def energy_single(i,j,k,Spin): #用于计算单个格点的能量 面心立方晶格的配位数为12
    theta1=Spin[i,j,k].real  #自旋的theta角
    phi1=Spin[i,j,k].imag  #自旋的phi角
    neighbour=[]  #用于存储相邻格点的位置
    #第一个面上相邻格点
    if j==N-1:
        neighbour.append((i,0,k))
    else:
        neighbour.append((i,j+1,k))
    if j==0:
        neighbour.append((i,N-1,k))
    else:
        neighbour.append((i,j-1,k))

    if i==N-1 and k!=0:
        neighbour.append((0,j,k-1))
    elif i !=N-1 and k==0:
        neighbour.append((i+1,j,N-1))
    elif i==N-1 and k==0:
        neighbour.append((0,j,N-1))
    else:
        neighbour.append((i+1,j,k-1))

    if i==0 and k!=N-1:
        neighbour.append((N-1,j,k+1))
    elif i !=0 and k==N-1:
        neighbour.append((i-1,j,0))
    elif i==0 and k==N-1:
        neighbour.append((N-1,j,0))
    else:
        neighbour.append((i-1,j,k+1))

    #第二个面上相邻格点
    if i==0:
        neighbour.append((N-1,j,k))
    else:
        neighbour.append((i-1,j,k))
    if i==N-1:
        neighbour.append((0,j,k))
    else:
        neighbour.append((i+1,j,k))
    
    if j==N-1 and k!=0:
        neighbour.append((i,0,k-1))
    elif j!=N-1 and k==0:
        neighbour.append((i,j+1,N-1))
    elif j==N-1 and k==0:
        neighbour.append((i,0,N-1))
    else:
        neighbour.append((i,j+1,k-1))

    if j==0 and k!=N-1:
        neighbour.append((i,N-1,k+1))
    elif j!=0 and k==N-1:
        neighbour.append((i,j-1,0))
    elif j==0 and k==N-1:
        neighbour.append((i,N-1,0))
    else:
        neighbour.append((i,j-1,k+1))
    #第三个面上相邻格点
    if k==0:
        neighbour.append((i,j,N-1))
    else:
        neighbour.append((i,j,k-1))
    if k== N-1:
        neighbour.append((i,j,0))
    else:
        neighbour.append((i,j,k+1))
    if i==N-1 and j!=0:
        neighbour.append((0,j-1,k))
    elif i!=N-1 and j==0:
        neighbour.append((i+1,N-1,k))
    elif i==N-1 and j==0:
        neighbour.append((0,N-1,k))
    else:
        neighbour.append((i+1,j-1,k))
    
    if i==0 and j!=N-1:
        neighbour.append((N-1,j+1,k))
    elif i!=0 and j==N-1:
        neighbour.append((i-1,0,k))
    elif i==0 and j==N-1:
        neighbour.append((N-1,0,k))
    else:
        neighbour.append((i-1,j+1,k))

    single_energy=0
    for a,b,c in neighbour:
        theta2=Spin[a,b,c].real
        phi2=Spin[a,b,c].imag
        single_energy += -cos(theta1)*cos(theta2)-sin(theta1)*sin(theta2)*(cos(phi1)*cos(phi2)+sin(phi1)*sin(phi2))
    return single_energy

def energy_total(Spin):  #用于计算体系整体的能量
    total_energy=0
    for i in range(N):
        for j in range(N):
            for k in range(N):
                total_energy+=energy_single(i,j,k,Spin)
    total_energy/=2  #每个相互作用能计算了两次，要除以2
    return total_energy

def random_spin(): #用于生成一个随机取向的自旋矢量
    phi=2*math.pi*random.uniform(0,1)
    theta=acos(1-2*random.uniform(0,1))
    return theta,phi

def Metropolis_MC(Spin,temperature):  #进行metropolis运算，随机改变一个格点自旋，选择性接受该调整
    i=random.randint(0,N-1)  #随机挑选出一个格点
    j=random.randint(0,N-1)
    k=random.randint(0,N-1)
    spin1=Spin[i,j,k]  #将原格点的状态存储
    energy1=energy_single(i,j,k,Spin)  #计算原格点的能量
    theta2,phi2=random_spin()  #生成一个新自旋矢量
    Spin[i,j,k]=theta2+1j*phi2
    energy2=energy_single(i,j,k,Spin)  #计算新格点的能量
    alpha=min(1,math.exp(-(energy2-energy1)/temperature))
    if random.uniform(0,1)>alpha:  #当概率满足这一条件是不接受，否则接受
        Spin[i,j,k]=spin1

def loop(N,temperature):  #对同一温度的点进行循环，使得能量达到稳定值
    Spin=np.zeros((N,N,N))*(1+0j)
    average_energy=0
    for i in range(max_loop_times):
        Metropolis_MC(Spin,temperature)
        if i > max_loop_times-100:
            average_energy+=energy_total(Spin)
    average_energy/=99
    return average_energy

for temperature in axi_tem:
    tem=loop(N,temperature)
    energy.append(tem)

plt.plot(axi_tem,energy)
plt.xlabel('T')
plt.ylabel('Energy')
plt.show()