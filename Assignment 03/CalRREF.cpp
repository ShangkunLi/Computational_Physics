/*
  Write a general code to transform a n*m matrix into the REDUCED ROW ECHELON FORM
  and use the code to obtain the RREF of the following matrix.
*/

/*
  本代码的基本实现思路为：
  1.首先将第1列的各元素a_i1(0<i<N+1)进行比较，将a_i1最大值所在的一行交换到第1行
  2.将第1行的元素全部除以a_11，将a_11的值变为1
  3.将第1行乘以某数之后加到后面各行，使得第1列中a_11以下各元素全部为0
  4.类似于上述步骤，对第2行的元素进行操作，使得第2行中a_22以下各元素全部为0
  5.依次，对第3、4...N行进行处理，将整个矩阵化为上三角矩阵
  6.再从最后一行开始，反向进行消元操作，将整个矩阵化为最简行阶梯形矩阵
*/

#include <iostream>

using namespace std;

int main()
{
    int i, j, k; //用于遍历矩阵
    int n,m; //用于存储矩阵的行数和列数
    double tmp;  //用于中间处理
    //由用户输入矩阵的行数n和列数m
    printf("\nPlease input the rows n and cols m of the Matrix:\n");
    scanf("%d %d",&n,&m);
    //创建一个n*m的矩阵Matrix[n][m];
    double** Matrix= new double* [n];
    for(i=0;i<n;i++)
    {
        Matrix[i]= new double [m];
    } 
    //按顺序输入矩阵的元素
    printf("\nPlease input %d numbers to form a %d*%d Matrix:\n", n * m, n, m);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            scanf(" %lf", &Matrix[i][j]);
        }
    }
    //打印初始矩阵
    printf("\nThe initial Matrix is:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%lf\t", Matrix[i][j]);
        }
        printf("\n");
    }
    //开始进行相关运算
    //从第一行开始，对每一行进行处理，将矩阵变为上三角矩阵
    for (i = 0; i < n; i++)
    {
        //首先将a[i][i]与a[j][i]进行比较，将max(a[j][i])的元素所在的行跟第i行交换位置
        for (j = i; j < n; j++)
        {
            if (Matrix[j][i] > Matrix[i][i])
            {
                swap(Matrix[j], Matrix[i]);
            }
        }
        tmp = Matrix[i][i];
        //进行操作，使矩阵变为上三角矩阵
        for (j = i; j < m; j++)
        {
            Matrix[i][j] /= tmp; //操作一，将第i行除以元素a[i][i]，使得元素a[i][i]为1
        }
        for (j = i + 1; j < n; j++)  //遍历第i行以下的各行，参数为j
        {
            tmp = Matrix[j][i];
            for (k = i; k < m; k++)  //遍历第j行的各个元素，参数为k
            {
                Matrix[j][k] -= Matrix[i][k] * tmp;  //操作二，用第i行的a[i][i]将a[j][i]消为0，i<j<N
            }
        }
    }
    //从下往上回代，将矩阵变为最简阶梯形矩阵
    for (i = n - 1; i >= 0; i--)
    {
        for (j = i - 1; j >= 0; j--)
        {
            tmp = Matrix[j][i];
            for (k = i; k < m; k++)
            {
                Matrix[j][k] -= Matrix[i][k] * tmp;
            }
        }
    }
    //打印出最终结果
    printf("\nThe reduced row echelon form is:\n");
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            printf("%lf\t", Matrix[i][j]);
        }
        printf("\n");
    }
    //删去给数组分配的空间
    for(i=0;i<n;i++)
    {
        delete[] Matrix[i];
    }
    delete[] Matrix;
    system("pause");
    return 0;
}