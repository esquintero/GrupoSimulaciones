#include "iostream"
#include "cmath"
using namespace std;

// Constante de tiempo de recuperacion fijo
const double Gamma = 0.08;
const double Tf = 7000;
const int N = 500; // Numero de betas
const double Betamin = 0.35; 
const double Betamax = 0.999; // Beta*s0 > Gamma
// Declaracion de funciones
double f1(double s, double i, double t);
double f2(double s, double i, double t);
void RK4(double & dt, double & s, double & i, double & t, double beta);

int main(void){
  cout.precision(8);
  // Condiciones iniciales
  double s = 0.999;
  double i = 0.001;
  double t = 0;
  double dt = 0.1;
  double Beta, dBeta; // Ahora varia Beta

  // Comportamiento de S_inf vs Ro
  double aux, R0;
  dBeta = (Betamax-Betamin)/N;
  for (int j = 0; j < N; ++j) {
    Beta = Betamin + j*dBeta;
    R0 = Beta/Gamma;
    while (t < Tf) {
      aux = s;
      RK4(dt,s,i,t,Beta);
      if (aux == s) {
	break;
      }
    }
    cout << R0 << " " << s << " " << Beta <<endl;
    t = 0;
    s = 0.999;
  }
    return 0;
}
  
// Implementando funciones
double f1(double s, double i, double t, double Beta)
{
  return -Beta*s*i;
}

double f2(double s, double i, double t, double Beta)
{
  return Beta*s*i - Gamma*i;
}

void RK4(double & dt, double & s, double & i, double & t, double beta)
{
  // Incrementos k para s y l para i
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;

  k1 = dt*f1(s,i,t,beta);
  l1 = dt*f2(s,i,t,beta);
  
  k2 = dt*f1(s + k1/2, i + l1/2, t + dt/2, beta);
  l2 = dt*f2(s + k1/2, i + l1/2, t + dt/2, beta);

  k3 = dt*f1(s + k2/2, i + l2/2, t + dt/2, beta);
  l3 = dt*f2(s + k2/2, i + l2/2, t + dt/2, beta);

  k4 = dt*f1(s + k3, i + l3, t + dt, beta);
  l4 = dt*f2(s + k3, i + l3, t + dt, beta);

  // Actualizar variables

  s += (k1 + 2*k2 + 2*k3 + k4)/6;
  i += (l1 + 2*l2 + 2*l3 + l4)/6;
  t += dt;

}

