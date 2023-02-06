/*
  Search for the minimum of the function
  g(x,y)=sin(x+y)+cos(x+2*y)
  in the whole place
*/
/*
  本题采用梯度下降法求得函数的全局最小值
  1.g(x,y)为一个周期函数，周期为2PI
  2.因此我们只需要在一个2PI*2PI的平面区域内找到最小值点，就可将其作为全局最小值点
  3.给定一个初始点(x0,y0),我们在[x0,x0+2PI]*[y0,yo+2PI]的平面内随机选取30个点
  4.从这30个点出发，分别沿负梯度方向移动，分别到达对应的极小值点
  5.在这30个局部极小值点中找到最小值，将其作为全局最小值点
  6.本算法中隐含一个假设：就是30次随机取点就一定可以找到函数所有的局部极小值
  该假设的合理性在于：平面的大小仅为2Pi*2PI，并且函数是连续的，因此选取30个点可以很好地代表这个局部平面内的点
*/
#include <iostream>
#include <cmath>
#include<ctime>
#include<cstdlib>

#define EPS 0.000001   //判断梯度为零时设置的精度值
#define PI 3.1415926

double grad[2]; //存储函数的梯度值

double g(double x, double y) //定义函数
{
  return sin(x + y) + cos(x + 2 * y);
}

void gra_g(double x, double y) //定义函数的梯度
{
  grad[0] = cos(x + y) - sin(x + 2 * y);
  grad[1] = cos(x + y) - 2 * sin(x + 2 * y);
}

//找到函数的全局最小值
void Find_min(double a, double x0, double y0)
{
  int i;
  double x, y;
  double min_x, min_y, minima;
  double min[30][2];
  srand(time(0));  //生成随机数种子
  for (i = 0; i < 30; i++)  //在2PI*2PI的平面上随机生成30个点，进行30次循环
  {
    x = x0 + 2.0 * PI * ((double)rand() / RAND_MAX);
    y = y0 + 2.0 * PI * ((double)rand() / RAND_MAX);
    gra_g(x, y);
    while (grad[0] > EPS || grad[1] > EPS)  //循环结束的标准，当梯度的两个分量都趋于0时，就认为到达了极小值点
    {
      x = x - a * grad[0];  //沿函数的负梯度方向移动
      y = y - a * grad[1];
      gra_g(x, y);
    }
    min[i][0] = x;  //将找到的函数极小值点放入一个2*2的二维数组中
    min[i][1] = y;
  }
  min_x = min[0][0];
  min_y = min[0][1];
  for (i = 1; i < 30; i++)  //从所有找到的极小值点中挑选出函数值最小的那一个
  {
    if (g(min[i][0], min[i][1]) < g(min_x, min_y))
    {
      min_x = min[i][0];
      min_y = min[i][1];
    }
  }
  minima = g(min_x, min_y);
  printf("The minima point is (%lf,%lf)\nThe minima is %lf\n", min_x, min_y, minima); //输出这个最小值
}

int main()
{
  double a, x0, y0;
  printf("Please input the step a: ");  //输入沿负梯度方向下降的步长，步长越小，精度越高
  scanf("%lf", &a);
  printf("Please input the initial location x0:\n");  //输入初始的点(x0,y0)，以该点为基准产生2PI*2PI的平面
  scanf("%lf", &x0);
  printf("Please input the initial location y0:\n");
  scanf("%lf", &y0);
  Find_min(a, x0, y0);
  system("pause");
  return 0;
}