#include <iostream>
#include <math.h>

using namespace std;

double f(double x)
{
    return (2 * x) - 5;
}

double gradient_descent(double x, double a)
{
    double tmp = x - (a * f(x));
    x = tmp;
    cout << x << endl;
    
    if (f(x) == 0)
        return x;
    else
        return gradient_descent(x, a);
}

int main()
{
    double learning_rate = 0.1;
    cout <<  gradient_descent(0, learning_rate) << endl;
}
