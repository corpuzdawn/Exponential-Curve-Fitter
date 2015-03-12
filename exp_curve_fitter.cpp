#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include "exp_curve_fitter.h"

using namespace std;

void expCurveFitter::importData(string filename )
{
  fstream f;
  double valueX;
  double valueY;
  f.open( filename.c_str() );

  int i=0;
  while (! f.eof() )
  {
    if ( f >> valueX )
      if (f >> valueY )
      {
        x_data.push_back( valueX);
        y_data.push_back( valueY);
      }
  }

}
/// NOTE: relative difference toleranc set at 0.01
void expCurveFitter::fit()
{
    int CurvatureStatus = this->determineCurve();
    ///initial A, B, C
    paramA = this->determineA( x_data, y_data);
    paramB = this->determineB(paramA, x_data, y_data);
    paramC = this->determineC(paramA, x_data, y_data);
    e0 = this->squareDiff(paramA, paramB, paramC);
  ///LOOP :: Calculate for  A0, A+, A- points for Aopt limited by relative difference
    while (this->relativeDifference() > 0.001){
    e0 = this->squareDiff(paramA, paramB, paramC);
    EParabolaX = this->EPointX(); ///still incomplete
    EParabolaY = this->EPointY(EParabolaX);
    paramA = this->EVertexA();
    this->showEData();
    paramB = this->determineB(paramA, x_data, y_data);
    paramC = this->determineC(paramA, x_data, y_data);
    e1 = this -> EVertexE();
    this->showE();
    }
    cout << "Final square difference: " << e0 <<endl;
}

int  expCurveFitter::determineCurve(){
    if ((y_data.at(1) - y_data.at(0)) > 0){
        if (y_data.at(2)-(2*y_data.at(1))+y_data.at(0) > 0) {
            return 1;
        }
        else
            return 4;
    }
    else {
        if (y_data.at(2)-(2*y_data.at(1))+y_data.at(0) > 0) {
            return 3;
        }
        else
            return 2;
            }
}

double expCurveFitter::determineA(vector<double> x_point, vector<double> y_point){
    double y=0;
    y=this->DetermineLowPoint(x_point, y_point);
    return (y-0.1*y);
}

double expCurveFitter::determineB(double A, vector<double> x_point, vector<double> y_point){
    double xy;
    double x;
    double y;
    double x2;
    vector<double> y_mod = y_point;

    ///Subtract A from Y values for new Y set
    for (int i = 0; i< y_point.size(); i++)
    {
        y_mod[i] = y_point.at(i)-A;
    }

    ///summation and linearization
    for (int i = 0; i< x_point.size(); i++)
    {
        xy += x_point.at(i) * log10(y_mod.at(i));
        x += x_point.at(i);
        y += log10(y_mod.at(i));
        x2 += x_point.at(i) * x_point.at(i);
    }

        return pow(10, (((y*x2)-(x*xy))/ ((x_point.size()*x2)-(x*x))));
}

double expCurveFitter::determineC(double A, vector<double> x_point, vector<double> y_point){
    double xy;
    double x;
    double y;
    double x2;
    vector<double> y_mod = y_point;

    ///Subtract A from Y values for new Y set
    for (int i = 0; i< y_point.size(); i++)
    {
        y_mod[i] = y_point.at(i)-A;
    }
    ///summation and linearization
     for (int i = 0; i< x_point.size(); i++)
    {
        xy += x_point.at(i) * log10(y_mod.at(i));
        x += x_point.at(i);
        y += log10(y_mod.at(i));
        x2 += x_point.at(i) * x_point.at(i);
    }
    return log(pow(10, ((x_point.size()* xy)- (x*y))/((x_point.size()*x2)-(x*x))));
}

///automatically determines A0,e0  A+,e+  A-,e-
vector<double> expCurveFitter::EPointX(){
    vector<double> points;
    points.push_back(paramA);
    points.push_back(paramA-increment);
    points.push_back(paramA+increment);
    incrementer();
    return points;
}

vector<double> expCurveFitter::EPointY(vector<double> A){
    vector<double> points;
    double B, C,e;
    for (int i =0; i < A.size(); i++){
        B = this->determineB(A.at(i), x_data, y_data);
        C = this->determineC(A.at(i), x_data, y_data);
        e = this->squareDiff(A.at(i),B,C);
        points.push_back(e);
    }
    return points;
}

/// NOTE determine e0,e-, e+,
double expCurveFitter::EVertexA(){
    double x1 = EParabolaX.at(0);
    double x2 = EParabolaX.at(1);
    double x3 = EParabolaX.at(2);
    double y1 = EParabolaY.at(0);
    double y2 = EParabolaY.at(1);
    double y3 = EParabolaY.at(2);

    double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	double A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	double B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	double C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

    return	-B / (2*A);
}

double expCurveFitter::EVertexE(){
    double x1 = EParabolaX.at(0);
    double x2 = EParabolaX.at(1);
    double x3 = EParabolaX.at(2);
    double y1 = EParabolaY.at(0);
    double y2 = EParabolaY.at(1);
    double y3 = EParabolaY.at(2);

    double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	double A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	double B     = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	double C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;

    return C - B*B / (4*A);
}

double expCurveFitter::relativeDifference(){
        return e0 - e1;
}

void expCurveFitter::incrementer(){
        increment  /= 100;
}

void expCurveFitter::resetIncrement(){
        increment  = 0.001;
}

double expCurveFitter::DetermineLowPoint(vector<double> x_point, vector<double> y_point){
    int counter = 0;
    double y = y_point.at(0);
    while(counter < y_point.size()){
            if (y > y_point.at(counter))
                y = y_point.at(counter);
            counter++;
    }
    return y;
}

int expCurveFitter::numPoints()
{
  return x_data.size();
}

void expCurveFitter::showData()
{
  for (int i=0; i<x_data.size(); i++)
  {
    cout << x_data.at(i)  << "\t\t" << y_data.at(i)  << endl;
  }
}

void expCurveFitter::showEData()
{
  cout << "A and e Parabolic points:" <<endl;
  for (int i=0; i<EParabolaX.size(); i++)
  {
    cout << EParabolaX.at(i)  << "\t\t" << EParabolaY.at(i)  << endl;
  }
}

void expCurveFitter::showE()
{
    cout << "e0: " << e0  << "\t\t e1: " << e1  << endl;
    cout << "relative difference: " << e0 - e1 <<endl;
}



double expCurveFitter::squareDiff(double A, double B, double C)
{
  double meanSquareDiff=0;
  double diff;

  for (int i=0; i<x_data.size(); i++ )
  {
    diff = A  + B* exp( C* x_data.at(i) ) - y_data.at(i);

    meanSquareDiff += diff*diff;
  }
  meanSquareDiff /= x_data.size();
  return  meanSquareDiff;
}

double expCurveFitter::getA()
{
  return paramA;
}

double expCurveFitter::getB()
{
  return paramB;
}

double expCurveFitter::getC()
{
  return paramC;
}
