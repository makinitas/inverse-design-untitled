#include <iostream>
#include <cmath>
#include <math.h>

using namespace std;

double f(double x)
{
    return (4 * pow(x, 3)) - 3;
}

double gradient_descent(double* x, double a)
{
    double tmp = *x - (a * f(*x));
    *x = tmp;

    //cout << *x << endl;
    //cout << round(f(*x) * 100000000.0) / 100000000.0 << endl << endl;
    
    if (round(f(*x) * 10000000000.0) / 10000000000.0 == 0)
        return *x;
    else
        return gradient_descent(x, a);
}

int main()
{
    double learning_rate = 0.1;
    double* x;
    *x = 0;
    cout <<  gradient_descent(x, learning_rate) << endl;
    //cout << f(0.909) << endl;
}
