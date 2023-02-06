'''
The interior of a 𝑑-dimensional hypersphere of unit radius is defined by the condition 𝑥_1^2+𝑥_2^2+⋯+𝑥_𝑑^2≤1. 
Write a program that finds the volume of a hypersphere using a Monte Carlo method. 
Test your program for 𝑑=2 and 𝑑=3 and then calculate the volume for 𝑑=4 and 𝑑=5.
Compare your results with the exact results. 
'''
import random
import math
import numpy as np

print('Please input the number of points:')
N=input()
N=int(N)

def hit_and_miss(n,N): #利用hit&miss方法来求解，n为空间维数，N为取点个数
    number=0
    for i in range(N):  #进行N次取点
        r=np.zeros(n)
        for j in range(n):  #对于每一个点随机生成n个坐标
            r[j]=random.uniform(-1,1)**2
        if sum(r)<=1:  #如果点在球内，则接受该点
            number+=1
    return (2**n*number/N)  #返回球的体积

for n in range(2,6):  #分别计算2-5维的超球体积
    print('\n维数：',n,'\t取点次数：',N,'\t体积计算值：',hit_and_miss(n,N),'\t体积理论值：',math.pow(math.pi,n/2)/math.gamma(n/2+1))