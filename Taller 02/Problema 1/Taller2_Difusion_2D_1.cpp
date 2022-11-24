#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "Random64.h"
#include <sstream>
#include <fstream>
using namespace std;

//ofstream outfile;

const int Lx=256, Ly=256; 
/*Lista de valores
1---- p_0=0.25, p=0.25
2---- p_0=0.1, p=0.2
3---- p_0=0.2, p=0.05
4---- p_0=0.15, p=0.1*/
const double p=0.1, p_0=0.15;//p_0 es la probabilidad de quedarse quieto y p la probabilidad de girar 90 grados

const int Q=4;//Número de flechas


//----------Clase LatticeGas------------
class LatticeGas{
private:
  int Vx[Q], Vy[Q]; //Vectores de cooredenadas de las direcciones 
  int n[Lx][Ly][Q], nnew[Lx][Ly][Q];//n[ix][i]
public:
  LatticeGas(void);
  void Borrese(void);
  void Inicie(int N, double mu, double sigma, Crandom &ran64);
  void Show(void);
  void Shownew(void);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double rho(int ix, int iy){ return (n[ix][iy][0]+n[ix][iy][1]+n[ix][iy][2]+n[ix][iy][3]);};
  double Varianza(void);
  void GrafiqueRho(void);
};


LatticeGas::LatticeGas(void){
  
  Vx[0]=1; Vx[1]=0; Vx[2]=-1; Vx[3]=0; 
  Vy[0]=0; Vy[1]=1; Vy[2]=0;  Vy[3]=-1;
  for(int ix = 0; ix < Lx; ix++)
	  for(int iy = 0; iy < Ly; iy++)
	    for(int i = 0; i < Q; i++) n[ix][iy][i] = nnew[ix][iy][i] = 0;

}

void LatticeGas::Borrese(void){
  for (int ix=0; ix<Lx;ix++){
    for(int iy=0; iy<Ly; iy++){
      for(int i=0;i<Q;i++){
	n[ix][iy][i]=nnew[ix][iy][i]=0;
      }
    }
  }
}

void LatticeGas::Inicie(int N, double mu, double sigma, Crandom &ran64){
  int ix,iy,i;
  
  while (N>0){ 
    ix=(int) ran64.gauss(mu,sigma);//Escoger una celda al azar
    iy=(int) ran64.gauss(mu,sigma);
    if(ix<0) ix=0; if(ix>Lx-1) ix=Lx-1;//Corregir en los bordes si es necesario
    if(iy<0) iy=0; if(iy>Ly-1) iy=Ly-1;
    i=(int) Q*ran64.r(); //Escoger una dirección al azar, elige un número entre 0 y Q
    
    if (n[ix][iy][i]==0) //si ese sitio está vacío
      { n[ix][iy][i]=1; N--;}//pongo una bolita ahí y decrezco N
    	 }
}


void LatticeGas::Show(void){
  for(int i=0;i<Q;i++){ 
    for (int ix=0; ix<Lx;ix++){
      for (int iy=0; iy<Ly;iy++){
	cout<<n[ix][iy][i];
      }
    }
    cout<<endl;
  }
}

void LatticeGas::Shownew(void){
  for(int i=0;i<Q;i++){
    for (int ix=0; ix<Lx;ix++){
      for (int iy=0; iy<Ly;iy++){
	cout<<nnew[ix][iy][i];
	
      }
    }
    cout<<endl;
  }
}




void LatticeGas::Colisione(Crandom & ran64){
  
  for (int ix=0; ix<Lx; ix++){ 
    for(int iy=0; iy<Ly; iy++){
double PB=ran64.r();//Generar un número al azar entre 0 y 1
      if(PB<=p_0){
	nnew[ix][iy][0]=n[ix][iy][0];//Girar 0°
	nnew[ix][iy][1]=n[ix][iy][1];
	nnew[ix][iy][2]=n[ix][iy][2];
	nnew[ix][iy][3]=n[ix][iy][3];
	}
      
      else if(PB>(p_0) && PB<=(p_0+p)){
	
	nnew[ix][iy][0]=n[ix][iy][3];//Gira 90° a la derecha
	nnew[ix][iy][1]=n[ix][iy][0];
	nnew[ix][iy][2]=n[ix][iy][1];
	nnew[ix][iy][3]=n[ix][iy][2];
      }

      else if(PB>(p_0+p) && PB<=(2*p+p_0)){
	nnew[ix][iy][0]=n[ix][iy][1]; //90 grados a la izquierda
	nnew[ix][iy][1]=n[ix][iy][2];
	nnew[ix][iy][2]=n[ix][iy][3];
	nnew[ix][iy][3]=n[ix][iy][0];
      }

      
      else{
	
	nnew[ix][iy][0]=n[ix][iy][2];//Gira 90° a la izquierda	
  nnew[ix][iy][1]=n[ix][iy][3];
	nnew[ix][iy][2]=n[ix][iy][0];
	nnew[ix][iy][3]=n[ix][iy][1];
      }
      
    
               
    }
  }
}

//Moverse a las siguientes celdas
void LatticeGas::Adveccione(void){
   for (int ix=0;ix<Lx;ix++){
     for (int iy=0;iy<Ly;iy++){
       for(int i=0;i<Q;i++){
	 n[(ix+Vx[i]+Lx)%Lx][(iy+Vy[i]+Ly)%Ly][i]=nnew[ix][iy][i];//Condiciones de frontera periodicas 
      }
     }
   }
}



double LatticeGas::Varianza(void){
  int ix, iy; double N, R,Nx,Ny, Rprom, xprom,yprom,Sigma2x,Sigma2y,Sigma2,SigmaR;
  
  //Calcular N
  for (N=0,ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      //N+=(n[ix][iy][0]+n[ix][iy][1]+n[ix][iy][2]+n[ix][iy][3]);
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


//-----------Programa Principal------------

int main(void){
  
  LatticeGas Difusion;
  Crandom ran64(8);
  int N=2400;
  double mu=Lx/2, sigma=16;
  int t, tmax=350;

  Difusion.Borrese();
  Difusion.Inicie(N, mu, sigma, ran64);
  ofstream outfile;
  outfile.open("T2_a_4.dat");

 for(t=0; t<tmax;t++){
   outfile<<t<<" "<<Difusion.Varianza()<<endl;
   Difusion.Colisione(ran64);
   Difusion.Adveccione();
   //Difusion.Show();
   //Difusion.Shownew();
  }
 //Difusion.GrafiqueRho();
 outfile.close();
 return 0;
}
