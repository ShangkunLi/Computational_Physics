!求解实系数一元二次方程的解
!系数a,b,c为任意实数，从键盘输入
program Solve_eq
implicit none
    real a, b, c, delta !分别声明方程系数和判别式
    real q !用于保证计算的精度
    real x1, x2 !声明方程的实数根
    real real_part, imag_part !声明方程复数根的实部和虚步
    print* , "Please input the coefficient:"
    read* , a, b, c !输入方程的系数
    !通过分析系数的取值，将方程分为了4种情况
    !1.所有系数全为零，方程恒成立，任意复数均为方程的根
    if (a == 0 .and. b==0 .and. c==0) then
        print'("Every complex number is root")'
    !2.a,b为零，常数项系数c不为零，方程不存在，方程无根
    else if (a == 0 .and. b==0 .and. c /=0) then
        print'("This equation does not have root.")'
    !3.二次项系数a为零，一次项系数b不为零，方程为一元一次方程，有一个根
    else if(a == 0 .and. b /= 0)then
        x1 = -c/b
        print'("This equation has one root:", f10.5)', x1
    !4.二次项系数a不为零，方程为二元一次方程，根据判别式判断其根的情况
    else
        delta = b*b-4.0*a*c !计算判别式
        !判别式大于零，方程有两个不相等的实根
        if (delta > 0.0) then
            q = (-b+sqrt(b*b-4.0*a*c))/2.0
            x1 = q/a
            x2 = c/q
            print '("This equation has two roots:", f10.5, f10.5)', x1, x2
        !判别式等于零，方程有两个相等的实根
        else if (delta == 0.0) then
            x1 = (-b+sqrt(b*b-4.0*a*c))/(2.0*a)
            print'("This equation has a double root:", f10.5)', x1
        !判别式小于零，方程有两个复根
        else
            real_part = -b/(2.0*a)
            imag_part = sqrt(4.0*a*c-b*b)/2.0*a
            print'("This equation has two complex roots:")'
            print'(f10.5,"+",f10.5,"i")', real_part, imag_part
            print'(f10.5,"-",f10.5,"i")', real_part, imag_part
        end if
    end if
    pause !为防止直接运行.exe文件时cmd在运行结束后直接关闭，特地增加此命令
end program Solve_eq