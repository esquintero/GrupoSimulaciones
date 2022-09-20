// Simular el movimiento de N granos con gravedad, fricción y choques inelásticos
#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
using namespace std;

//---- declarar constantes ---
const double g=9.8, K=1.0e4, Gamma=50;
const double Lx=60, Ly=120;
const int Nx=5, Ny=5, N=Nx*Ny;

//Constantes de la Fuerza de Lennard Jones
const double E=1.0, r0=1.0;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//--- declarar clases -----
class Cuerpo;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau; double I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
	      double theta0,double omega0);
  void BorreFuerza(){F.load(0,0,0); tau=0;};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void AdicioneTorque(double tau0){tau+=tau0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
  double Gettheta(void){return theta;}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,double m0,double R0,
		    double theta0,double omega0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt);  theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m);  omega+=tau*(Coeficiente*dt/I);
}
/*void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
  // cout<<" , "<<r.x()<<"+"<<R*cos(theta)/7<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7<<"*t";
  }*/


//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Grano);
  void CalculeFuerzaEntre(Cuerpo & Grano1, Cuerpo & Grano2);
  void CalculeFuerzaPared(Cuerpo & Grano1);
};

void Colisionador::CalculeFuerzas(Cuerpo * Grano){
  int i,j; vector3D Fg;
  //--- Borrar todas las fuerzas ---
  for(i=0;i<N;i++)
    Grano[i].BorreFuerza(); 
  //--- Sumar el peso ---
  for(i=0;i<N;i++){
    Fg.load(0,-Grano[i].m*g,0); 
    Grano[i].AdicioneFuerza(Fg);
    CalculeFuerzaPared(Grano[i]);
  }
  //--- Calcular Fuerzas entre pares de granos ---
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++)
      CalculeFuerzaEntre(Grano[i], Grano[j]);
}


void Colisionador::CalculeFuerzaPared(Cuerpo &Grano1){
  double x=Grano1.Getx(), y=Grano1.Gety();
  vector3D r=Grano1.r;
  double h, K=1.0e4, d=r.norm(), R=Grano1.R;
  vector3D n;
  
    //Pared de la izquierda
  
      if(x<R){
    h=R-x;
    n.load(1,0,0);
    vector3D F=n*(K*pow(h,1.5));  
    Grano1.AdicioneFuerza(F);
      }

      //Pared de la derecha

      if((x+R)>Lx){
	h=R-(Lx-x);
    n.load(1,0,0);
    vector3D F=n*(K*pow(h,1.5));  
    Grano1.AdicioneFuerza(F*(-1));
      }

      //Pared de abajo

      if(y<R){
	 h=R-y;
    n.load(0,1,0);
    vector3D F=n*(K*pow(h,1.5));  
    Grano1.AdicioneFuerza(F);
      }


      //Pared de arriba

      if((y+R)>Ly){
	h=R-(Ly-y);
    n.load(0,1,0);
    vector3D F=n*(K*pow(h,1.5));  
    Grano1.AdicioneFuerza(F*(-1));
      }
      
}
  


//Fuerza entre moléculas de Lennard Jones

void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2){
  vector3D r21,n,F2; double d21,F;
  r21=Grano2.r-Grano1.r; d21=r21.norm(); n=r21/d21;
  F=12*E/d21*(pow(r0/d21,12)-pow(r0/d21,6));
  F2=F*n; Grano2.AdicioneFuerza(F2); Grano1.AdicioneFuerza(F2*(-1));
}





//----------------- Funciones de Animacion ----------
/*void InicieAnimacion(void){
  // cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'Gas2D.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}*/

//-----------  Programa Principal --------------  
int main(void){
  Cuerpo Grano[N];
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=2.5, kT=10, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,tmax=200,/*tmax=10*(Lx/V0),*/ tcuadro=tmax/1000,dt=0.001;
  //double dx=Lx/(Nx+1), dy=Ly/(Ny+1);
  double Theta;
  double y, yprom;
  
  //InicieAnimacion(); //Dibujar
  

  /*//Inicializar las paredes
  
  double Rpared=100*Lx, Mpared=100*m0; //radio y masa de las paredes
  //---------------(  x0,       y0,Vx0,Vy0,    m0,    R0, theta0,omega0) 
  Grano[N+0].Inicie(Lx/2,Ly+Rpared,  0,  0,Mpared,Rpared,      0,     0); //Pared de arriba
  Grano[N+1].Inicie(Lx/2,  -Rpared,  0,  0,Mpared,Rpared,      0,     0); //Pared de abajo
  Grano[N+2].Inicie(Lx+Rpared,Ly/2,  0,  0,Mpared,Rpared,      0,     0); //Pared derecha
  Grano[N+3].Inicie(  -Rpared,Ly/2,  0,  0,Mpared,Rpared,      0,     0); //Pared izquierda
  */
  
  
  //Inicializar las moléculas
  
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();//el ángulo de cada molécula respecto a x es aleatorio 
      //--------------------(   x0,   y0,          Vx0,          Vy0, m0,R0,theta0,omega0)
      Grano[Nx*iy+ix].Inicie((ix+1)*10,(iy+1)*10,V0*cos(Theta),V0*sin(Theta), m0,R0,0,1);
    }
  
  for(t=0/*tdibujo=0*/ ; t<tmax ; t+=dt/*,tdibujo+=dt*/){
    //Dibujar
    y=0;
    yprom=0;
    //if(tdibujo>tcuadro){
     
      //InicieCuadro();
      for(i=0;i<N;i++) {//Grano[i].Dibujese();
	y+=Grano[i].Gety();
	}
      yprom=y/N;
      cout<<t<<"\t"<<yprom<<"\n";
      
      //TermineCuadro();
      
      //tdibujo=0;
      //}

    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Grano);
    for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Grano[i].Mueva_r(dt,epsilon);  

  }   

  
  
  return 0;
}

  
