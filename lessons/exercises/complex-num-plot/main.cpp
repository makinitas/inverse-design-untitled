#include <iostream>
#include <math.h>
#include <string>
#include "../../../utils/pbPlots/pbPlots.h"
#include "../../../utils/pbPlots/supportLib.h"

#define PI 3.14159265

using namespace std;

void output_data(double x[], double y[], int n)
{
    bool success;
    StringReference *errorMessage = new StringReference();
	RGBABitmapImageReference *imageReference = CreateRGBABitmapImageReference();

    vector<double> xs;
	vector<double> ys;

    for (int i = 0; i < (2 * n) + 1; i++)
    {
        xs.push_back(x[i]);
        ys.push_back(y[i]);
    }

	success = DrawScatterPlot(imageReference, 600, 400, &xs, &ys, errorMessage);

    if(success){
        vector<double> *pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, "plots/plot" + to_string(n) + ".png");
        DeleteImage(imageReference->image);
	}
}

double calculate_complex_num(int n, double ang)
{
    double fr, fi;

    for (int l = -n; l <= n; l++)
    {
        fr += cos(l*ang);
        fi += sin(l*ang);
    }

    return ((fr * fr) + (fi * fi)) / pow((2 * n) + 1, 2);
}

int main(int argc, char** argv)
{
    for (int j = 3; j < 100; j++)
    {
        int n = j;
        double y[(2 * n) + 1];
        double x[(2 * n) + 1];
    
        // populate arrays
        for (int i = -n; i <= n; i++)
        {
            x[n + i] = i * 0.01;
            y[n + i] = calculate_complex_num(n, x[n + i]);
            cout << x[n + i] << "   " << y[n + i] << endl;
        }
    
        // plot graph
        output_data(x, y, n);
    }
    
    return 0;
}
