/*
  Proyecto Final:
  Lattice Boltzmann para la ecuación AD en coordenadas cartesianas con una fuente
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
using namespace std;

ofstream outfile;

//-----Constantes Globales-----//
const int Lx = 200;                 //tamaño de la simulacion
const int Ly = 200;

const int Q = 9;                    //numero de direcciones

const double C = 1.0/sqrt(3.0);   //Velocidad C correcta para D2Q9
//const double C = 1.01;            //Velocidad necesaria para que la varianza se comporte como difusión
const double C2 = C*C;

const double Dcoef = 0.05;
const double tau = Dcoef/C2 + 0.5;
const double Utau = 1.0/tau;
const double UmUtau = 1.0-Utau;

//-----Clase LatticeBoltzmann------
class LatticeBoltzmann{
private:
  double w[Q];        //pesos por direccion
  int Vx[Q], Vy[Q];   //vectores de velocidad
  double *f, *fnew;   //funciones de distribucion - * es un apuntador que lo lleva a la direccion de memoria
  double Snew[Lx][Ly][Q], Sold[Lx][Ly][Q];
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};         //indice de las celdas
  double rho(int ix, int iy, bool UseNew);
  double Ssum(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double S(int ix, int iy, int t);
  double Si(int ix, int iy, double Ux0, double Uy0, int t, int i);           //forzamiento
  double feq(double rho0, double Ux0, double Uy0, int i);
  void Start(double rho0, double Ux0, double Uy0);
  void Collision(double Ux0, double Uy0);
  void ImposeFields(void);
  void Advection(double Ux0, double Uy0, int t);    //Advección Camilo
  double Varianza(void);
  void Print(const char * NameFile);
  /*double rekt(void){
    int i; double sum;
    for(sum=0,i=0;i<Q;i++){
      sum+=w[i]*Vy[i]*Vy[i];
    }
    return sum;
  };*/ //Función para rectificar que los pesos y velocidades estén bien, devolviendo 1/3
};
//Constructor
LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
 w[0]=4.0/9.0;  w[1]=w[2]=w[3]=w[4]=1.0/9.0;  w[5]=w[6]=w[7]=w[8]=1.0/36.0;
 //Set the velocity vectors
 Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
 Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;

           Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
           Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;
  //Crea arrays dinamicos
  int ArraySize = Lx*Ly*Q;
  f = new double [ArraySize];  fnew = new double [ArraySize];
}
//Destructor
LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
}
//Rho
double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}
//Ssum
double LatticeBoltzmann::Ssum(int ix, int iy, bool UseNew){
  double sum; int i;
  for(sum=0,i=0;i<Q;i++){
    if(UseNew) sum+=Snew[ix][iy][i]; else sum+=Sold[ix][iy][i];
  }
  return sum;
}
//Componentes del vector J
double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0,i=0;i<Q;i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0,i=0;i<Q;i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}
//Trve Forzamiento
double LatticeBoltzmann::S(int ix, int iy, int t){
  if(iy==100 && ix==100){
    /*double dens=0;
    if(t<=5){ dens=1.0; }  //Fuente que alimenta durante 10 pasos de tiempo y luego se apaga
    return dens;*/
    return 1.0;
  }
  //Autopistas en el sur
  /*if(iy==120 && (ix >0 && ix<132)){
      return 1;
  }
  if((ix==45 || ix ==90) && (iy>0 && iy<120)){
      return 1;
  }
  if(iy==60 && (ix >45 && ix<90)){
      return 1;
  }
  //kennedy<<
  if(ix==60 && (iy >120 && iy<270)){
      return 1;
  }
  if(ix==20 && (iy >120 && iy<200)){
      return 1;
  }
  if(ix==110 && (iy >120 && iy<180)){
      return 1;
  }
  if(ix==150 && (iy >120 && iy<170)){
      return 1;
  }
  if(iy==200 && (ix >20 && ix<60)){
      return 1;
  }
  if(iy==170 && (ix >60 && ix<150)){
      return 1;
  }
  //teusaquillo, engativá y chapinero
  if(ix==80 && (iy >150 && iy<310)){
      return 1;
  }
  if(ix==130 && (iy >150 && iy<270)){
      return 1;
  }
  if(iy==180 && (ix >80 && ix<130)){
      return 1;
  }
    if(iy==220 && (ix >0 && ix<130)){
      return 1;
  }
    if(iy==250 && (ix >70 && ix<170)){
      return 1;
  }
  if(iy==310 && (ix >70 && ix<180)){
      return 1;
  }
  if(iy==295 && (ix >100 && ix<170)){
      return 1;
  }
    if(ix==80 && (iy >220 && iy<310)){
      return 1;
  }
  if(ix==95 && (iy >220 && iy<310)){
      return 1;
  }
  //autonorte
  if(ix==145 && (iy >220 && iy<400)){
      return 1;
  }
  if(ix==175 && (iy >220 && iy<400)){
      return 1;
  }
  if(ix==100 && (iy >220 && iy<350)){
      return 1;
  }
  //Fabricas
  if(ix==138 && iy==270){
    return 1;            // Bavaria cra 53 # 127
  }
  if(iy==177 && ix==50){
    return 1;            // Industria Nacional de Gaseosas cll 25 # 95
  }
  if(iy==100 && ix==80){
    return 1;            // General motors cll 56sur # 36
  }
  if(iy==150 && ix==160){
    return 1;            // Diana cra 13 #93
  }
  if(iy==245 && ix==145){
    return 1;            // Nestle diag 92 # (cra)19
  }*/
  else {
    return 0;
  }
}
//Forzamiento LBGK
double LatticeBoltzmann::Si(int ix, int iy, double Ux0, double Uy0, int t, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i];
  return w[i]*S(ix,iy,t)*(1+(UdotVi/C2));
}
//Funcion equilibrio
double  LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  //FEQUILIBRIO TESIS JULIANA #1
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1.0+(UdotVi/C2)+(pow(UdotVi,2)/(2*pow(C2,2)))-(U2/(2*C2)));
}
//Start
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix, iy, i, n0; //double Sold, Snew;
  for(ix=0; ix<Lx; ix++)    //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){   //en cada direccion
        n0 = n(ix,iy,i);
        f[n0] = feq(rho0,Ux0,Uy0,i);
        //Calculamos S con t=0 y t=1
        //Debido a que S es indep. de t tendrán los mismos valores todo el tiempo
        Sold[ix][iy][i] = Si(ix,iy,Ux0,Uy0,0,i);
        Snew[ix][iy][i] = Si(ix,iy,Ux0,Uy0,1,i);
      }
}
//Colision
void LatticeBoltzmann::Collision(double Ux0, double Uy0){
  int ix,iy,i,n0; double rho0;//,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++){
      //Calcule los campos macroscopicos en la celda
      rho0 = rho(ix,iy,false);
      //Ux0 = 0; Uy0 = 0;
      for(i=0;i<Q;i++){     //para cada vector de velocidad
        n0 = n(ix,iy,i);
        fnew[n0] = UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i)+(3.0*Snew[ix][iy][i])/2.0-Sold[ix][iy][i]/2.0;
      }
    }
}
//Imponer campos para cerrar las fronteras
void LatticeBoltzmann::ImposeFields(void){
  int ix,iy,i,n0;
  for(ix=0,iy=0;iy<Ly;iy++) //Pared izquierda
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
  for(ix=Lx-1,iy=0;iy<Ly;iy++) //Pared derecha
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
  for(ix=0,iy=Ly-1;ix<Lx;ix++) //Pared arriba
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
  for(ix=0,iy=0;ix<Lx;ix++) //Pared abajo
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
}
//ADVECCION ORIGINAL DE CAMILO
void LatticeBoltzmann::Advection(double Ux0, double Uy0, int t){
  int ix, iy, i, ixnext, iynext, ixback, iyback, n0, n0next;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){     //en cada direccion
        if( (ix==0) || (ix==Lx-1) || (iy==0) || (iy==Ly-1) ){ continue; } //No hay adveccion en fronteras
        ixnext = (ix+Vx[i]+Lx)%Lx; iynext = (iy+Vy[i]+Ly)%Ly;
        //faltaría un if que corte las paredes (?)
        ixback = (ix-Vx[i]+Lx)%Lx; iyback = (iy-Vy[i]+Ly)%Ly;
        n0 = n(ix,iy,i); n0next = n(ixnext,iynext,i);
        f[n0next] = fnew[n0];     //fronteras periodicas
        //Sold = Snew;
        //Parece que cambiar los valores de Snew y S no tiene efecto porque siempre se mantienen iguales
        Sold[ix][iy][i] = Snew[ix][iy][i];
        Snew[ix][iy][i] = Si(ix,iy,Ux0,Uy0,t+1,i); //Debe ser t+1 porque este Snew será el del siguiente ciclo
        //Sold[ix][iy][i] = Snew[ixback][iyback][i];
        //Snew[ix][iy][i] = Si(ix,iy,Ux0,Uy0,t,i);
      }
}
double LatticeBoltzmann::Varianza(void){
  int ix, iy; double N,Sigma2x,Sigma2y,SigmaR,xprom,yprom;

  //Calcular Ntotal
  for(N=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      N+=rho(ix,iy,true);
    }
  }
  N=N-(Lx*Ly);
  //Calcular yprom
  for(yprom=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      yprom+=iy*(rho(ix,iy,true)-1.0); //Este -1 ya asegura que el prom se mide perfectamente
    }
  }
  yprom/=N;
  //calc xprom
  for(xprom=0,iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
      xprom+=ix*(rho(ix,iy,true)-1.0); //Este -1 ya asegura que el prom se mide perfectamente
    }
  }
  xprom/=N;

  //Calcular Sigma2
  for(Sigma2y=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      Sigma2y+=pow((iy-yprom),2)*(rho(ix,iy,true)-1.0);
    }
  }
  Sigma2y/=(N-1);

  for(Sigma2x=0,iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
      Sigma2x+=pow((ix-xprom),2)*(rho(ix,iy,true)-1.0);
    }
  }
  Sigma2x/=(N-1);

for( SigmaR=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      SigmaR+=(pow((iy-yprom),2)+(pow((ix-xprom),2.0)))*(rho(ix,iy,true)-1.0);
    }
  }
  SigmaR/=(N-1);

  return SigmaR;
}
//Print
void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix, iy;
  for(ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      rho0=rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//-----Programa Principal-----

int main(void){
  LatticeBoltzmann Ondas;
  int t, tmax = 200;
  double rho0 = 1.0, Ux0 = 0.0, Uy0 = 0.0;
  outfile.open("varianza.dat");

  Ondas.Start(rho0, Ux0, Uy0);
  //Evolucione
  for(t=1; t<tmax; t++){
    if(t>=50){outfile<<t<<" "<<Ondas.Varianza()<<endl;} //Escribimos desde t=50 para la varianza
    Ondas.Collision(Ux0,Uy0);
    Ondas.ImposeFields();
    Ondas.Advection(Ux0, Uy0, t);
  }
  outfile.close();
  //Imprima
  Ondas.Print("D2Q9.dat");

  return 0;
}
