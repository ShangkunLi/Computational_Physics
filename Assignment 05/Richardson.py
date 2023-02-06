'''
Compute the derivative of f (x) = sin x at x = π/3 
using the Richardson extrapolation algorithm. 
Start with h = 1 and find the number of rows in the Richardson table 
required to estimate the derivative with six significant decimal digits. 
Output the Richardson table.
'''
'''
本题目利用Richardson方法进行一阶导数的计算
1. 定义Richardson方法中需要的D(n,0)函数:D_n0()
2. 定义一个矩阵MatrixD，用于存储Richardson table
3. 从D(0,0)出发，这相当于一个1*1矩阵
4. 先往下添加一个元素D(1,0),再由D(1,0)和D(0,0)计算出D(1，1),这就得到了一个2*2矩阵
5. 再往下添加一个元素D(2,0)，重复上述步骤，得到D(2,2)
6. 依次类推，每算出一个对角元素D(i,i)时，将其与D(i-1,i-1)比较一下，若精度满足要求，则输出结果
'''
import math

h=1.0  #设置h的初始值
pi=math.pi
eps=1e-6  #设置精度
MatrixD=[]  #声明一个矩阵，用于存储Richardson table

def f(x):  #定义待求导函数
    return math.sin(x)

def D_n0(n):  #定义D(n,0)
    return (f(pi/3+h/2**n)-f(pi/3-h/2**n))/(h/2**(n-1))

i=1
MatrixD.append([D_n0(0)])
while True:
    MatrixD.append([])  #新增一行
    MatrixD[i].append(D_n0(i))  #计算下一行的D(n,0)
    for j in range(1,i+1):    #计算下一行的D(n,1),D(n,2),...,D(n,n)
        MatrixD[i].append(MatrixD[i][j-1]+(MatrixD[i][j-1]-MatrixD[i-1][j-1])/(4**j-1))
    if abs(MatrixD[i][i]-MatrixD[i-1][i-1])<eps:  #当满足精度要求时，退出循环
        break
    i+=1  #每经历一次循环，i的序号加1

print('The number of rows in the Richardson table required to \
estimate the derivative with six significant decimal digits is',len(MatrixD),'\n')
print('The derivative of f(x) at x = π/3 is %.6f \n' %MatrixD[-1][-1])
print('The Richardson table is \n')
for i in range(len(MatrixD)):
    print(MatrixD[i])