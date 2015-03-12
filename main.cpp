#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "exp_curve_fitter.h"

using namespace std;

int main(int argc, char *argv[])
{
    string s;
    expCurveFitter fitter;

    cout << "*************************************" <<endl
            << "******Exponential Curve Fitter*******" <<endl
            << "***CE 30**Corpuz**Del Rosario**Lim***" <<endl
            << "*************************************" <<endl
            << "\nThis program accepts a set of (x,y) points"
            <<"\nand accurately fits the data with an exponential curve. " <<endl;


    cout << "\nfilename to import data from:";
    getline(cin, s);
    cout << "Data set: " <<endl;
    cout << " x\t\t y" << endl;
    fitter.importData(s);
    fitter.showData();
    cout << "Curvature type: " <<fitter.determineCurve()  <<endl;
    fitter.fit();
    cout <<endl;
    cout << "Exponential Fit Equation: " << fitter.getA() << " + ";
    cout << fitter.getB() << " ";
    cout << "exp ( " << fitter.getC() << " x) " <<endl;
    cout << "Press the enter key to continue ...";
    cin.get();
    return EXIT_SUCCESS;
}
