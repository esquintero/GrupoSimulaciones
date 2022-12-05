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
const int Lx = 100;                 //tamaño de la simulacion
const int Ly = 100;

const int Q = 9;                    //numero de direcciones
const double W0 = 1.0/3.0;          //cte que define los pesos

const double C = 0.6;               //velocidad de la onda ajustada a 0.6- por estabilidad numerica: C < 0.707 cells/click
const double C2 = C*C;

const double tau = 0.8;             //Valor para que D=0.05
const double Utau = 1.0/tau;
const double UmUtau = 1-Utau;

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
  void Collision(void);
  //void ImposeFields(int t, double Ux0, double Uy0);
  //void ImposeFields(void);
  void Advection(double Ux0, double Uy0, int t);    //Advección Camilo
  //void Advection(void);                               //Advección Profe original
  double Varianza(void);
  void Print(const char * NameFile);
};
//Constructor
LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
 w[0]=4.0/9;  w[1]=w[2]=w[3]=w[4]=1.0/9;  w[5]=w[6]=w[7]=w[8]=1.0/36;
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
  double sum; int i,n0;
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
  if(iy==50 && ix==50){
      return 2;
  }
  else {
    return 0;
  }
  /*//Autopistas en el sur
  if(iy==120 && (ix >0 && ix<132)){
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
  }   //Set the weights
  w[0]=4.0/9;  w[1]=w[2]=w[3]=w[4]=1.0/9;  w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;

            Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
            Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;

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
      return 1;Fx
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
  }1
  if(iy==150 && ix==160){
    return 1;            // Diana cra 13 #93
  }
  if(iy==245 && ix==145){
    return 1;            // Nestle diag 92 # (cra)19
  }*/
}
//Forzamiento LBGK
double LatticeBoltzmann::Si(int ix, int iy, double Ux0, double Uy0, int t, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i];
  return w[i]*S(ix,iy,t)*(1+(UdotVi/C2));
}
//Funcion equilibrio
double  LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  //FEQUILIBRIO TESIS JULIANA
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1.0+(UdotVi/C2)+(pow(UdotVi,2)/(2*pow(C2,2)))-(U2/(2*C2)));
  //FEQUILIBRIO D2Q9 DADA POR PROFE EN CLASE
  //double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  //return rho0*w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
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

  /*//PONEMOS UN RHO DISTINTO AL INICIO EN LA CELDA CENTRAL
  //Calcule los campos macroscopicos en la celda
  double rhoC = 3.0;
  for(i=0;i<Q;i++){     //para cada vector de velocidad
    n0 = n(50,50,i);
    f[n0] = feq(rhoC,Ux0,Uy0,i);
    fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
  }*/
}
//Colision
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++){
      //Calcule los campos macroscopicos en la celda
      rho0 = rho(ix,iy,false);
      Ux0 = Jx(ix,iy,false)/rho0; Uy0 = Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){     //para cada vector de velocidad
        n0 = n(ix,iy,i);
        fnew[n0] = UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i)+(3.0*Snew[ix][iy][i])/2.0-Sold[ix][iy][i]/2.0;
      }
    }
}
//Imponer campos
/*/void LatticeBoltzmann::ImposeFields(int t, double Ux0, double Uy0){
void LatticeBoltzmann::ImposeFields(){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  ix=50,iy=50;
  //Calcule los campos macroscopicos en la celda
  rho0 = 3.0; Ux0=0; Uy0=0;
  for(i=0;i<Q;i++){     //para cada vector de velocidad
    n0 = n(ix,iy,i);
    f[n0] = feq(rho0,Ux0,Uy0,i);
    fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
  }
}*/
//ADVECCION ORIGINAL DE CAMILO
void LatticeBoltzmann::Advection(double Ux0, double Uy0, int t){
  int ix, iy, i, ixnext, iynext, ixback, iyback, n0, n0next;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){     //en cada direccion
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
/*/ADVECCION ORIGINAL D2Q9 DADA EN CLASE PROFE
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}*/
double LatticeBoltzmann::Varianza(void){
  int ix, iy; double N,R,Rprom,Sigma2x,Sigma2y,SigmaR,xprom,yprom;

  //Calcular Ntotal
  /*for(N=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      N+=rho(ix,iy,true);
    }
  }*/
  N=2.0;
  //Calcular yprom
  for(yprom=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      yprom+=iy*rho(ix,iy,true);
    }
  }
  yprom/=N;
  //calc xprom
  for(xprom=0,iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
      xprom+=ix*rho(ix,iy,true);
    }
  }
  xprom/=N;

  //Calcular Sigma2
  for(Sigma2y=0,ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      Sigma2y+=(pow((iy-yprom),2.0)*rho(ix,iy,true));
    }
  }
  Sigma2y/=(N-1);

  for(Sigma2x=0,iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
      Sigma2x+=(pow((ix-xprom),2.0)*rho(ix,iy,true));
    }
  }
  Sigma2x/=(N-1);

for( SigmaR=0, ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      SigmaR+=((pow((iy-yprom),2.0)+(pow((ix-xprom),2.0)))*rho(ix,iy,false));
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
  int t, tmax = 20;
  double rho0 = 1.0, Ux0 = 0.0, Uy0 = 0.0;
  //outfile.open("dftcont.dat");

  Ondas.Start(rho0, Ux0, Uy0);
  //Evolucione
  for(t=1; t<tmax; t++){
    //outfile<<t<<" "<<Ondas.Varianza()<<endl;
    Ondas.Collision();
    //Ondas.ImposeFields();
    Ondas.Advection(Ux0, Uy0, t); //Adveccion de Camilo
    //Ondas.Advection();          //Advección Profe original
  }
  //outfile.close();
  //Imprima
  Ondas.Print("D2Q9.dat");

  return 0;
}
