#include "iostream"
#include "cmath"
using namespace std;

// Constantes de contagio y tiempo de recuperacion
const double Beta = 0.35;
const double Gamma = 0.08;
// const double Gamma = 0.36; //Punto b
const double Tf = 70;

// Declaracion de funciones
double f1(double s, double i, double t);
double f2(double s, double i, double t);
double f3(double i, double t);
void RK4(double & dt, double & s, double & i, double & t);
double Recuperados(double & t, double & i, double & dt);

int main(void){
  // Condiciones iniciales
  double s = 0.999;
  double i = 0.001;
  double r = 0; // Al principio no hay recuperados
  double dt = 0.01;

  // Epidemia
  for ( double t = 0; t < Tf; ){
    cout << t << " " << s << " " << i << " " << r << " " << s + i + r << endl;
    RK4(dt, s, i, t);
    r += Recuperados(t, i, dt);
  }
  
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

double f3(double i, double t)
{
  return Gamma*i;
}

void RK4(double & dt, double & s, double & i, double & t)
{
  // Incrementos k para s y l para i
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;

  k1 = dt*f1(s,i,t);
  l1 = dt*f2(s,i,t);
  
  k2 = dt*f1(s + k1/2, i + l1/2, t + dt/2);
  l2 = dt*f2(s + k1/2, i + l1/2, t + dt/2);

  k3 = dt*f1(s + k2/2, i + l2/2, t + dt/2);
  l3 = dt*f2(s + k2/2, i + l2/2, t + dt/2);

  k4 = dt*f1(s + k3, i + l3, t + dt);
  l4 = dt*f2(s + k3, i + l3, t + dt);

  // Actualizar variables y calculo de r

  s += (k1 + 2*k2 + 2*k3 + k4)/6;
  i += (l1 + 2*l2 + 2*l3 + l4)/6;
  t += dt;
}

double Recuperados(double & t, double & i, double & dt)
{
  double k1,k2,k3,k4;
  
  k1 = dt*f3(i,t);
  k2 = dt*f3(i,t+dt/2);
  k3 = dt*f3(i,t+dt/2);
  k4 = dt*f3(i,t+dt);

  return (k1 + 2*k2 + 2*k3 + k4)/6;
}
