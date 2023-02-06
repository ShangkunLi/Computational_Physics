'''
The table below gives the temperature T along a metal rod whose ends are kept at fixed constant temperature.
The temperature is a function of the distance x along the rod.
(1)Compute a least-squares, straight-line fit to these data using T(x)=a+bx
(2)Compute a least-squares, parabolic-line fit to these data using T(x)=a+bx+cx^2
'''
'''
本题目利用最小二乘法进行求解
通过分析可知满足ka^2最小时，系数应当是一个线性方程组的解
通过求该线性方程组即可求出系数的取值
'''

import numpy as np
import matplotlib.pyplot as plt

#输入原始数据
x=[1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0]
T=[14.6,18.5,36.6,30.8,59.2,60.1,62.2,79.4,99.9]
#声明方程Ax=b中的系数矩阵A和向量b
#Matrix1和Vector1分别对应线性拟合时的系数矩阵和向量b
Matrix1=np.zeros((2,2))
Vector1=[0.0,0.0]
#Matrix2和Vector2分别对应二次拟合时的系数矩阵和向量b
Matrix2=np.zeros((3,3))
Vector2=[0.0,0.0,0.0]
#通过已有数据给系数矩阵和向量赋值
for i in range(9):
    Matrix1[0][0]+=1
    Matrix1[0][1]+=x[i]
    Matrix1[1][0]+=x[i]
    Matrix1[1][1]+=x[i]**2
    Vector1[0]+=T[i]
    Vector1[1]+=T[i]*x[i]

    Matrix2[0][0]+=1
    Matrix2[0][1]+=x[i]
    Matrix2[0][2]+=x[i]**2
    Matrix2[1][0]+=x[i]
    Matrix2[1][1]+=x[i]**2
    Matrix2[1][2]+=x[i]**3
    Matrix2[2][0]+=x[i]**2
    Matrix2[2][1]+=x[i]**3
    Matrix2[2][2]+=x[i]**4
    Vector2[0]+=T[i]
    Vector2[1]+=T[i]*x[i]
    Vector2[2]+=T[i]*x[i]**2
#求解线性方程
Solution1=np.matmul(np.linalg.inv(Matrix1),Vector1)
a1=Solution1[0]
b1=Solution1[1]
Solution2=np.matmul(np.linalg.inv(Matrix2),Vector2)
a2=Solution2[0]
b2=Solution2[1]
c2=Solution2[2]
#打印结果
print('Straight-line Fit: T=a+bx  ','a=',a1,'b=',b1)
print('Parabolic-line Fit: T=a+bx+cx^2  ','a=',a2,'b=',b2,'c=',c2)
#绘制出拟合得到的函数
x_set=np.linspace(0,10,200)
T1_set=[a1+b1*x for x in x_set]
T2_set=[a2+b2*x+c2*x**2 for x in x_set]
plt.scatter(x,T,color='g')
plt.plot(x_set,T1_set,color='r',label='Straight-line Fit')
plt.plot(x_set,T2_set,color='b',label='Parabolic-line Fit')
plt.xlabel('x/cm')
plt.ylabel('T/C')
plt.legend()
plt.show()

