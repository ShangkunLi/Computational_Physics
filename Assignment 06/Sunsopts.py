'''
Download the file called sunspots.txt, which contains the observed number of sunspots on the Sun for each month since January 1749.
Write a program to calculate the Fourier transform of the sunspot data
and then make a graph of the magnitude squared |ck|^2 of the Fourier coefficients
as a function of k—also called the power spectrum of the sunspot signal. 
You should see that there is a noticeable peak in the power spectrum at a nonzero value of k. 
Find the approximate value of k to which the peak corresponds.
What is the period of the sine wave with this value of k?
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft

#声明两个数组，用于存储sunspots.txt中的数据
x_list=[]
y_list=[]
#读取文件，并将数据存储到声明的两个数组中
file_obj=open("sunspots.txt",mode='r')
for line in file_obj.readlines():
    line=line.strip()
    x,y=line.split('\t')
    x_list.append(x)
    y_list.append(y)

N=len(x_list)
y_fft=abs(fft(y_list))  #进行快速傅立叶变换

#变换之后对数据进行归一化，求出幅值；进行取半操作，因为FFT的结果关于频率中心对称，故还需进行取半操作
y_fft_normalization=y_fft/(N/2)  #归一化操作
y_fft_normalization[0] /= 2
x_list_half = x_list[0:int(N/2)]  #取半操作
y_fft_normalization_half=y_fft_normalization[0:int(N/2)]

plt.xticks(np.arange(-1,3000,500))
plt.xlabel('k')
plt.ylabel('|Ck^2|')
plt.plot(x_list_half[1:],(y_fft_normalization_half[1:])**2)  #作出能谱图，其中去掉了无效的零频数据
plt.show()
plt.xlabel('k')
plt.ylabel('|Ck^2|')
plt.xticks(np.arange(-1,50,10))
plt.plot(x_list_half[1:50],(y_fft_normalization_half[1:50])**2)
plt.show()
