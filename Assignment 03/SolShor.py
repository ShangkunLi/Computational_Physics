'''
  Solve the 1D Schrodinger equation with the potential
  (i) V(x)= x^2; (ii) V(x)= x^4-x^2
  with the variational approach using a Gaussian basis (either fixed widths or fixed centers).
  Consider the three lowest energy eigenstates.
'''
'''
  整体的思路为：
  1.固定高斯函数的v=0.5，改变s的值，由此计算出对应的哈密顿矩阵H和矩阵S
  2.求解广义本征值问题HC=ESC
  3.对求得的本征值E进行排序，得到最小的三个即为对应的能量值
'''
import numpy as np
import math
import scipy.linalg as linalg
import matplotlib.pyplot as plt

#输入用于求解的高斯基的个数
print("Plese input the number of Gaussian basis:")
n = input()
n = int(n)

#定义高斯基函数，取定v=0.5，参数为s和x
def Gauss(x,si):
  return (math.sqrt((0.5/np.pi)))*math.exp(-0.5*(x-si)**2)

#对于V(x)=x^2，定义其哈密顿矩阵H1的第i行j列元素的计算公式
def H1_ij(si, sj):
    return (1.0+si*sj) * math.exp(-0.25 * (si - sj)**2) / (2.0 * math.sqrt(np.pi))

#对于V(x)=x^4-x^2，定义其哈密顿矩阵H2的第i行j列元素的计算公式
def H2_ij(si, sj):
    return (12.0+pow(si,4)+4.0*pow(si,3)*sj+4.0*(sj**2)+pow(sj,4)+4.0*si*sj*(6.0+sj**2)+(si**2)*(4.0+6.0*(sj**2))) * math.exp(-0.25 * (si - sj)**2) / (32.0 * math.sqrt(math.pi))

#对于矩阵S，定义其第i行j列元素的计算公式
def S_ij(si, sj):
    return math.exp(-0.25 * (si - sj)**2) / (2.0*math.sqrt(math.pi))

#定义一个绘图函数，用于绘制对应势能的最低三个能量的解
def Plot_eigenstate(sort_eigen,eigenstate):
  psidata1=[]  #用于存储波函数的取值
  psidata2=[]
  psidata3=[]
  axis=np.linspace(-4,4,200)
  for x in axis:
    psi1=0
    psi2=0
    psi3=0
    for i in range(n):  #依次计算三个能量本征值的波函数
      si = -((n-1)/4)+0.5*i
      psi1 += eigenstate[i][sort_eigen[0]]*Gauss(x,si)
      psi2 += eigenstate[i][sort_eigen[1]]*Gauss(x,si)
      psi3 += eigenstate[i][sort_eigen[2]]*Gauss(x,si)
    psidata1.append(psi1)
    psidata2.append(psi2)
    psidata3.append(psi3)
  plt.plot(axis,psidata1)  #绘制图像
  plt.plot(axis,psidata2)
  plt.plot(axis,psidata3)
  plt.xlabel('x')
  plt.ylabel('psi')

#定义矩阵的维数
H1 = np.zeros((n, n), dtype=float)
H2 = np.zeros((n, n), dtype=float)
S = np.zeros((n, n), dtype=float)

#生成n*n的哈密顿矩阵H和矩阵S，并赋值
for i in range(n):
    si = -((n-1)/4)+0.5*i
    for j in range(n):
      sj = -((n-1)/4)+0.5*j
      H1[i,j] = H1_ij(si, sj)
      H2[i,j] = H2_ij(si, sj)
      S[i,j] = S_ij(si, sj)

#求解广义本征值问题HC=ESC，E用于存储本征值，eigenstate用于存储本征向量
E1,eigenstate1=linalg.eig(H1,S)
E2,eigenstate2=linalg.eig(H2,S)
E1=abs(E1)
E2=abs(E2)
#对能量进行从小到大的排序
sort_eigen1=E1.argsort()
sort_eigen2=E2.argsort()

print("\nThe three lowest energy eigenstates of potential V(x)=x^2 are:")
#对其能量值进行排序，取最小三个能量值输出
print(E1[sort_eigen1][0:3])
print("\nThe three lowest energy eigenstates of potential V(x)=x^4-x^2 are:")
#对其能量值进行排序，取最小三个能量值输出
print(E2[sort_eigen2[0:3]])
print('\n')
#输出势能V(x)=x^2与V(x)=x^4-x^2时的波函数图像
Plot_eigenstate(sort_eigen1,eigenstate1)
plt.title('The three lowest wave function of potential x^2')
plt.show()
Plot_eigenstate(sort_eigen2,eigenstate2)
plt.title('The three lowest wave function of potential x^4-x^2')
plt.show()