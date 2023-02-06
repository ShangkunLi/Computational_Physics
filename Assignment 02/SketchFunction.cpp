/*
  Sketch the function x^3-5x+3=0
  (1) Determine the two positive roots to 4 decimal places using the bisection
  method. Note: You first need to bracket each of the roots.
*/
/*本过程由7个函数组成
  1.f函数定义一元三次方程对应的函数
  2.f_1函数为f函数的导数
  3.FindBracket函数用于找到两个正根分别所在的小区间（这一步也可用观察函数图像代替，但C/C++中进行函数画图不方便，所以采取此方法）
  4.FindRoot_bis函数利用二分法找到方程在每个小区间中的解
  5.FindRoot_new函数用牛顿迭代法找到方程靠近初值的解
  6.FindRoot_hyb函数利用混合方法找到方程的解
  7.main函数执行所有功能
*/

#include <iostream>
#include <cmath>
#include<cstdlib>

// 1.定义需要求解方程对应的函数f
double f(double x) { return x * x * x - 5 * x + 3; }
// 2.定义函数f的导函数f_1
double f_1(double x) { return 3 * x * x - 5; }

// 3.定义函数FindBracket，用于找到函数两个正根所在的大致区间
//这一步也可用观察函数图像代替，但C/C++中进行函数画图不方便，所以采取此方法
void FindBracket(double bracket[][2])
{
    int i = 0;
    float start = 0, end = 10;    //粗略估计方程的两个正根都位于0到10之间
    for (start = 0; start < end;) //因为已知函数有两个正根，因此以0.5为步长，找出使得函数符号相反的两个区间
    {
        if (f(start) * f(start + 0.5) < 0)
        {
            bracket[i][0] = start; //将区间存储在二维数组bracket中，其中数组每一行为一个区间
            bracket[i][1] = start + 0.5;
            i += 1;
        }
        start += 0.5;
    }
}

// 4.定义函数FindRoot_bis，用二分法找方程的根，输入值为根所在的大致区间
float FindRoot_bis(double a, double b)
{
    while (fabs(a - b) > 0.00001) //精度控制，当区间长度小于0.00001时跳出循环，这就保证了结果精确到四位小数
    {
        if (f(a) * f((a + b) / 2) < 0)
        {
            b = (a + b) / 2; //如果a端函数值与中点函数值异号，就将b的值取为中点值
        }
        else
        {
            a = (a + b) / 2; //如果a端函数值与中点函数值同好，就将a的值取为中点值
        }
    }
    return a; //当精度满足要求时，返回所求的根
}

// 5.定义函数FindRoot_new，使用牛顿迭代法求根，输入值为根附近的一个值
double FindRoot_new(double x)
{
    while (fabs(f(x) / f_1(x)) > pow(10, -14)) //控制精度，当对x的修正值小于10^(-14)时，这就保证了结果精确到十四位小数
    {
        x = x - f(x) / f_1(x);
    }
    return x; //当精度满足要求时，返回所求的根
}

// 6.定义函数FindRoot_hyb,使用混合方法求根，输入值为根所在的大致区间
double FindRoot_hyb(double a, double b)
{
    double x;
    x = (a + b) / 2;                  //先取a，b的中点作为第一个根
    while (fabs(f(x)) > pow(10, -14)) //控制精度，当某个根所对应的函数值小于10^(-14)时，就认为精度满足要求
    {
        if (f_1(x) != 0) //如果切线并不平行于x轴，就采用牛顿迭代法进行下一次迭代
        {
            x = x - f(x) / f_1(x);
            if (x > a && x < b) //用牛顿迭代法所求得的根如果在当前区间内，就根据这个根的函数值取新的区间
            {
                if (f(a) * f(x) < 0)
                {
                    b = x; //如果a点函数值与x点函数值符号相反，就将x取为b
                }
                else
                {
                    a = x; //如果a点函数值与x点函数值符号相同，就将x取为a
                }
            }
            else
            {
                x = (a + b) / 2; //用牛顿迭代法所求得的根如果不在当前区间内，舍弃这个根，然后采用二分法求解
                if (f(a) * f(x) < 0)
                {
                    b = x;
                }
                else
                {
                    a = x;
                }
                x = (a + b) / 2;
            }
        }
        else //如果切线平行于x轴，就采用二分法进行下一次迭代
        {
            x = (a + b) / 2;
            if (f(a) * f(x) < 0)
            {
                b = x;
            }
            else
            {
                a = x;
            }
            x = (a + b) / 2;
        }
    }
    return x; //当精度满足要求时，返回所求的根
}

int main()
{
    double x1, x2;
    double bracket[2][2];                            //用于存储两个正根分别所在的小区间
    FindBracket(bracket);                            //找到这两个正根分别所在的小区间
    x1 = FindRoot_bis(bracket[0][0], bracket[0][1]); //采用二分法求解两个正根
    x2 = FindRoot_bis(bracket[1][0], bracket[1][1]);
    printf("\nUsing the bisection method:\nx1=%.4f\nx2=%.4f\n", x1, x2); //以四位小数输出到屏幕
    x1 = FindRoot_new(x1);                                               //采用牛顿迭代法求解两个正根
    x2 = FindRoot_new(x2);
    printf("\nUsing the Newton-Rapson method:\nx1=%.14f\nx2=%.14f\n", x1, x2); //以使十四位小数输出到屏幕
    x1 = FindRoot_hyb(bracket[0][0], bracket[0][1]);                           //采用混合法求解两个正根
    x2 = FindRoot_hyb(bracket[1][0], bracket[1][1]);
    printf("\nUsing the Hybrid method:\nx1=%.14f\nx2=%.14f\n\n", x1, x2); //以十四位小数输出到屏幕
    system("pause");
    return 0;
}