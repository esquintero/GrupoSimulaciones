#include <iostream>
#include <cmath>
#include "../vector.h"
#include "planetas_lib.h"
#include <fstream>
/////////////////////////////////////////////////////////////////
// Cuerpos Celestes - metodo de integracion PEFRL              //
//                                                             //
// En este programa se simula y anima el movimiento ente dos   //
// planetas con un observador ajeno al sistema, a este último  //
// se le llama colisionador o interaccionador.                 //
// Los objetos celestes son el Sol y Júpiter a quienes vamos a //
// asignar las masas y radios como razones entre ellos:        //
//                                                             //
// m_sol = 1047 m_jupiter     ,  r_sol = 9.73 r_jupiter        //
// m_jupiter = 1              ,  r_jupiter = 1                 //
//                                                             //
/////////////////////////////////////////////////////////////////
const double G=1.0;
// Para el análisis de Sol-Júpiter hacemos N=2.
const int N=2;

// Coeficientes de PEFRL
const double Zeta   = 0.1786178958448091e00;
const double Lambda = -0.2123418310626054e00;
const double Chi    = -0.6626458266981849e-1;
const double Coeficiente1 =(1-2*Lambda)/2;
const double Coeficiente2 =1-2*(Chi+Zeta);


int main(){
  Cuerpo Planeta[N];
  Colisionador Newton;
  double m1=1047, m2=1, r=1000;
  double theta=M_PI/3;
  double M=m1+m2, x1=-m2*r/M, x2=m1*r/M;
  double omega = sqrt(G*M/pow(r,3.0)), T=2*M_PI/omega, V1=omega*x1, V2=omega*x2;
  double t, dt=0.01;
  double tmax=21.0*T;
  double tdibujo, tcuadro=T/100;
  int i;
  std::ofstream outfile;

  outfile.open("S-J_original_axis.dat");

  Planeta[0].Inicie(x1, 0,0,  0, V1,0,m1,97.3); // Sol
  Planeta[1].Inicie(x2, 0,0,  0, V2,0,m2,10.0); // Júpiter
  for(t=0,tdibujo=0; t<tmax; t+=dt, tdibujo+=dt){
    outfile<<Planeta[0].Getx()<<" "<<Planeta[0].Gety()<<" "
	     <<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<" "
	     <<Planeta[2].Getx()<<" "<<Planeta[2].Gety()<<std::endl;
    // Haga el movimiento y los calculos por PEFRL
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Zeta);            // 1 
    Newton.CalculeFuerza(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);     // 2
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Chi);             // 3
    Newton.CalculeFuerza(Planeta);     
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);           // 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Coeficiente2);     // 5
    Newton.CalculeFuerza(Planeta);     
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Lambda);           // 4
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Chi);             // 3
    Newton.CalculeFuerza(Planeta);
    for(i=0;i<N;i++) Planeta[i].Mueva_v(dt,Coeficiente1);     // 2
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt, Zeta);            // 1
  }
  return 0;
}
