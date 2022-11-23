#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <sstream>
#include <fstream>
#include "Random64.h"

using namespace std;

ofstream outfile;


const int Lx=256, Ly=256;
const double p0=0.25, p=0.25;
//p_0 es la probabilidad de quedarse quieto y p la probabilidad de girar a la derecha 90 grados

const int Q=4;//Número de flechas


//----------Clase LatticeGas------------
class LatticeGas{
private:
  int Vx[Q], Vy[Q]; //Vectores de cooredenadas de las direcciones 
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q];
public:
  LatticeGas(void);
  void Inicie(int N, double mu, double sigma);
  void Show(void);
  void Shownew(void);
  void Colisione(void);
  void Adveccione(void);
  double rho(int ix, int iy); 
  double Varianza(void);
  void GrafiqueRho(void);
};


LatticeGas::LatticeGas(void){
  Vx[0]=1; Vx[1]=0; Vx[2]=-1; Vx[3]=0; 
  Vy[0]=0; Vy[1]=1; Vy[2]=0;  Vy[3]=-1;
}


void LatticeGas::Inicie(int N, double mu, double sigma){
  for(int ix=0;ix<Lx;ix++){
    for(int iy=0;iy<Ly;iy++){
      double rho=(N/(sigma*sigma*2.0*M_PI))*(exp(-0.5*pow((ix-mu)/sigma,2.0)-0.5*pow((iy-mu)/sigma,2.0)));// rho es la densidad de probabilidad Gaussiana
      for(int i=0;i<Q;i++)
      f[ix][iy][i]=rho/Q;
    }
  }
}



void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){ 
    for (int ix=0; ix<Lx;ix++){
      for (int iy=0; iy<Ly;iy++){
	outfile<<f[ix][iy][i];
      }
    }
    outfile<<endl;
  }
}

void LatticeGas::Shownew(void){
  for(int i=0;i<Q;i++){ 
    for (int ix=0; ix<Lx;ix++){
      for (int iy=0; iy<Ly;iy++){
	outfile<<fnew[ix][iy][i];
      }
    }
    outfile<<endl;
  }
}


void LatticeGas::Colisione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
	fnew[ix][iy][i]=p0*f[ix][iy][i]+p*f[ix][iy][(i+1)%4]+(1-p0-2*p)*f[ix][iy][(i+2)%4]+p*f[ix][iy][(i+3)%4];
}

//Moverse a las siguientes celdas
void LatticeGas::Adveccione(void){
  for (int ix=0; ix<Lx;ix++){
    for (int iy=0; iy<Ly; iy++){
      for(int i=0;i<Q;i++){
	f[(ix+Vx[i]+Lx)%Lx][(iy+Vy[i]+Ly)%Ly][i]=fnew[ix][iy][i];//Condiciones de frontera periodicas 
      }
    }
  }
}


double LatticeGas::rho(int ix, int iy){
  double suma; int i;
  for(suma=0, i=0; i<Q; i++)
    suma+=f[ix][iy][i];
  return suma;
}


double LatticeGas::Varianza(void){
  int ix, iy; double N,Nx,Ny, R, Rprom, Sigma2x,Sigma2y,SigmaR,xprom,yprom;
  
  //Calcular N
  for (N=0,ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      N+=rho(ix,iy);
    }
  }
  
  //Calcular yprom
  for(yprom=0,Ny=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      yprom+=iy*rho(ix,iy);
    }
  }
  yprom/=N;
  //calc xprom
  for(xprom=0,Nx=0, iy=0; iy<Ly; iy++){
    for(ix=0; ix<Lx; ix++){
      xprom+=ix*rho(ix,iy);
    }
  }
  xprom/=N;
  
  //Calcular Sigma2
  for(Sigma2y=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      Sigma2y+=(pow((iy-yprom),2.0)*rho(ix,iy));
    }
  }
  Sigma2y/=(N-1);

  for(Sigma2x=0, iy=0; iy<Ly; iy++){
    for(ix=0; ix<Lx; ix++){
      Sigma2x+=(pow((ix-xprom),2.0)*rho(ix,iy));
    }
  }
  Sigma2x/=(N-1);

for( SigmaR=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      SigmaR+=((pow((iy-yprom),2.0)+(pow((ix-xprom),2.0)))*rho(ix,iy));
    }
  }
  SigmaR/=(N-1);
  return SigmaR; 
}

void LatticeGas::GrafiqueRho(void){
  for(int ix=0;ix<Lx;ix++)
    for(int iy=0;iy<Ly;iy++)
      outfile<<ix<<iy<<" "<<rho(ix,iy)<<endl;
}


//-----------Programa Principal------------

int main(void){
  
  LatticeGas Difusion;
  int N=2400;
  double mu=Lx/2, sigma=16;
  int t, tmax=350;
  outfile.open("dftcont.dat");
  
  Difusion.Inicie(N, mu, sigma);
  
  for(t=0; t<tmax;t++){
    outfile<<t<<" "<<Difusion.Varianza()<<endl;
    Difusion.Colisione();
    Difusion.Adveccione();
    //Difusion.Shownew();
  }
 // Difusion.GrafiqueRho();
 // Difusion.Show();
 outfile.close();
 return 0;
}
