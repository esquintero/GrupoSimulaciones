#include "iostream"
#include "cmath"
using namespace std;

// Constantes de contagio y tiempo de recuperacion
const double Beta = 0.35;
const double Gamma = 0.08;

// Declaracion de funciones
double f1(double s, double i, double t);
double f2(double s, double i, double t);

int main(void){
  return 0;
}

// Implementando funciones
double f1(double s, double i, double t)
{
  return -Beta*s*i;
}

double f2(double s, double i, double t)
{
  return Beta*s*i - Gamma*i;
}

