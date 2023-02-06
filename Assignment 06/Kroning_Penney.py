'''
One-dimensional Kronig-Penney problem
Periodic potential V(x)=V(x+a)
Using FFT, find the lowest three eigenvalues of the eigenstates that satisfy:
pesai(x)=pesai(x+a)
U0=2eV, Lw=0.9nm, LB=0.1nm
'''
'''
本题的主要求解思路是：
1.先将势函数进行傅立叶分解
2.然后由分解得到的势函数的分立的幅值去计算哈密顿矩阵
3.最后求解哈密顿矩阵的本征值问题，得到本征值，即为对应的
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
import scipy.linalg as linalg

def V(x):  #定义周期性的势函数
    if 0<=x<=0.90:
        return 2
    return 0

N=2000 #定义取样点的个数

x_list=np.arange(0,1,1/N)  #将势函数的横坐标分为N份
V_list=[V(x) for x in x_list]  #取对应的势函数的值，即取样

V_fft_central=abs(fft(V_list))  #对势函数进行快速傅立叶变换，得出其系数的模
V_fft_central /= (N/2)  #对模进行归一化处理，得到对应的幅值
V_fft_central[0]/=2

#由于一scipy库中的FFT函数结果是关于中间点对称，因此我们需要进行一个变换，将FFT的结果变成关于原点对称
V_fft_origin=np.append(V_fft_central[int(N/2):],V_fft_central[0:int(N/2)])

H=np.zeros((N,N))  #构建哈密顿算符对应的矩阵
for p in range(N):  #公式推导中p的取值范围为-N/2,...,N/2，此处一一对应为0,1,...,N来代替
    for q in range(N):
        if 0<=p-q<=N-1:
            H[p][q]=V_fft_origin[p-q]
        if p==q:  #为简化数值计算，我们取h^2/(2*m*a^2)=1
            H[p][q]+=(q-N/2)**2  #此处减去N/2的原因是，公式推导中的指标与此处的求和指标相差N/2

eigen_val=abs(linalg.eig(H)[0])  #求解本征值问题
eigen_val.sort()   #对本征值进行排序
print('The three lowest eigenvalues are:')
print (eigen_val[1:7])   #因为能量存在简并，故输出前6个本征值，第1个本征值为0，无物理意义，故舍去
plt.xlabel('n')
plt.ylabel('En')
plt.plot(eigen_val[1:10])
plt.show()
