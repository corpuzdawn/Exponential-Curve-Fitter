#ifndef EXPCURVEFITTER
#define EXPCURVEFITTER 1
#include <fstream>
#include <string>
#include <vector>
#include "exp_curve_fitter.h"

using namespace std;

class expCurveFitter  // fits data to A + B exp ( C x)
{
  public:
    void importData(string filename );
    void fit();
    double squareDiff( double A, double B, double C);
    double getA();
    double getB();
    double getC();
    void showData();
    int numPoints();
    double getX(int index);
    double getY(int index);

    ///project methods
    int determineCurve();
    double determineA(vector<double> x_point, vector<double> y_point);
    double determineB(double A, vector<double> x_point, vector<double> y_point);
    double determineC(double A, vector<double> x_point, vector<double> y_point);
    double DetermineLowPoint(vector<double> x_point, vector<double> y_point);
    vector<double> EPointX();
    vector<double> EPointY(vector<double> A);
    double EVertexA();
    double EVertexE();
    double relativeDifference();
    void incrementer();
    void resetIncrement();
    void showEData();
    void showE();

  private:
    vector<double> x_data;
    vector<double> y_data;
    vector<double> EParabolaX;
    vector<double> EParabolaY;
    double paramA;
    double paramB;
    double paramC;
    double e0 = 0; ///initial value
    double e1 = 0; ///initial value
    double increment= 0.001; ///initial value

};
#endif
