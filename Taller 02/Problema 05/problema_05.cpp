#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#define Lx 1
#define Ly 64
#define N 32 //Threads per Block
const int M=(Lx*Ly+N-1)/N; //Blocks per Grid

#define Q 9

const int ArraySize=Lx*Ly*Q;

//Más adelante serán definidas en la matriz tau[4]
//Tau de 1/2 hace que colapse el LB
const float tau=1.2;       //tau[0]
const float Utau=1.0/tau;  //tau[1]
const float UmUtau=1-Utau; //tau[2]
const float ThreeUmU2tau=3*(1-1/(2*tau)); //tau[3]

//------------ PROGRAMMING ON THE DEVICE ----------------
//---------------Constants (Symbols)----------------
__constant__ float d_w[9];      //Weights
__constant__ int d_Vx[9];       //Velocity vectors
__constant__ int d_Vy[9];
__constant__ float d_tau[4];

//---------- DEVICE FUNCTIONS ----
__device__ int d_n(int ix,int iy,int i){
  return (ix*Ly+iy)*Q+i;
}
__device__ float d_rho(int ix,int iy,float *d_f,float *d_fnew){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=d_n(ix,iy,i);
    sum+=d_f[n0];
  }
  return sum;
}
__device__ float d_Jx(int ix,int iy,float Fx,float *d_f,float *d_fnew){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=d_n(ix,iy,i);
    sum+=d_Vx[i]*d_f[n0];
  }
  return sum+0.5*Fx;
}
__device__ float d_Jy(int ix,int iy,float Fy,float *d_f,float *d_fnew){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=d_n(ix,iy,i);
    sum+=d_Vy[i]*d_f[n0];
  }
  return sum+0.5*Fy;
}
__device__ float d_Fi(float Ux0,float Uy0,float Fx,float Fy,int i){
  float UdotVi=Ux0*d_Vx[i]+Uy0*d_Vy[i];
  float FdotVi=Fx*d_Vx[i]+Fy*d_Vy[i], UdotF=Ux0*Fx+Uy0*Fy;
  return d_tau[3]*d_w[i]*(FdotVi-UdotF+3*UdotVi*FdotVi);
}
__device__ float d_feq(float rho0,float Ux0,float Uy0,int i){
  float UdotVi=Ux0*d_Vx[i]+Uy0*d_Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*d_w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
//---------- KERNELS -------------
__global__ void d_Collision(float d_gx,float d_gy,float *d_f,float *d_fnew){
  int icell,ix,iy,i,n0;  float rho0,Ux0,Uy0; float Fx,Fy;
  //Find which thread an which cell should I work
  //ix=0;
  icell=blockIdx.x*blockDim.x+threadIdx.x;
  ix=icell/Ly; iy=icell%Ly;
  //compute the macroscopic fields on the cell
  rho0=d_rho(ix,iy,d_f,d_fnew); Fx=d_gx*rho0; Fy=d_gy*rho0;
  Ux0=d_Jx(ix,iy,Fx,d_f,d_fnew)/rho0; Uy0=d_Jy(ix,iy,Fy,d_f,d_fnew)/rho0;
  for(i=0;i<Q;i++){ //for each velocity vector
    n0=d_n(ix,iy,i);
    d_fnew[n0]=d_tau[2]*d_f[n0]+d_tau[1]*d_feq(rho0,Ux0,Uy0,i)+d_Fi(Ux0,Uy0,Fx,Fy,i);
    }
  }
/*__global__ void d_Collision(float d_gx,float d_gy,float *d_f,float *d_fnew,float *d_matoi){
    //Find which thread an which cell should I work
    //int ix,iy;
    //Find which thread an which cell should I work
    ix=0;
    iy=blockIdx.x*blockDim.x+threadIdx.x;
    d_matoi[iy]=1.0;
}*/
__global__ void d_ImposeFields(float *d_f,float *d_fnew,float rho_bot,float rho_top){
  int i,iy,n0;
  //Lower wall
  iy=0;
  for(i=0;i<Q;i++){n0=d_n(0,iy,i); d_fnew[n0]=d_feq(rho_bot,0,0,i);}
  //Upper wall
  iy=Ly-1;
  for(i=0;i<Q;i++){n0=d_n(0,iy,i); d_fnew[n0]=d_feq(rho_top,0,0,i);}
}
__global__ void d_Advection(float *d_f,float *d_fnew){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  ix=0;
  iy=blockIdx.x*blockDim.x+threadIdx.x;
  for(i=0;i<Q;i++){ //on each direction
    ixnext=(ix+d_Vx[i]+Lx)%Lx; iynext=(iy+d_Vy[i]+Ly)%Ly;
    n0=d_n(ix,iy,i); n0next=d_n(ixnext,iynext,i);
    d_f[n0next]=d_fnew[n0]; //periodic boundaries
  }
}
//-------------------- CLASE: LatticeBoltzmann -----------
class LatticeBoltzmann{
private:
  float h_w[Q];         //Weights
  int h_Vx[Q],h_Vy[Q];  //Velocity vectors
  float h_tau[4];         //Constants:
  float *h_f, *h_fnew; float *d_f, *d_fnew; //Distribution Functions
  //float *h_matoi, *d_matoi;
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int h_n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  float h_rho(int ix,int iy,bool UseNew);
  float h_Jx(int ix,int iy,bool UseNew,float Fx);
  float h_Jy(int ix,int iy,bool UseNew,float Fy);
  float h_Fi(float Ux0,float Uy0,float Fx,float Fy,int i);
  float h_feq(float rho0,float Ux0,float Uy0,int i);
  void Start(float rho0,float Ux0,float Uy0);
  void Collision(float gx,float gy);
  void ImposeFields(void);
  void Advection(void);
  void Print(const char * NameFile,float gx,float gy);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //CONSTANTS(d_Symbols)
  //---Charge constantes on the Host-----------------
  //Taus
  h_tau[0]=tau;  h_tau[1]=Utau;  h_tau[2]=UmUtau;  h_tau[3]=ThreeUmU2tau;
  //Set the weights
  h_w[0]=4.0/9;  h_w[1]=h_w[2]=h_w[3]=h_w[4]=1.0/9;  h_w[5]=h_w[6]=h_w[7]=h_w[8]=1.0/36;
  //Set the velocity vectors
  h_Vx[0]=0;  h_Vx[1]=1;  h_Vx[2]=0;  h_Vx[3]=-1; h_Vx[4]=0;
  h_Vy[0]=0;  h_Vy[1]=0;  h_Vy[2]=1;  h_Vy[3]=0;  h_Vy[4]=-1;

              h_Vx[5]=1;  h_Vx[6]=-1; h_Vx[7]=-1; h_Vx[8]=1;
              h_Vy[5]=1;  h_Vy[6]=1;  h_Vy[7]=-1; h_Vy[8]=-1;
  //------Send to the Device-----------------
  cudaMemcpyToSymbol(d_w,h_w,Q*sizeof(float),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vx,h_Vx,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_Vy,h_Vy,Q*sizeof(int),0,cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(d_tau,h_tau,4*sizeof(float),0,cudaMemcpyHostToDevice);
  //Create the dynamic arrays
  h_f=new float [ArraySize];  h_fnew=new float [ArraySize];
  //h_matoi=new float [64];
  //Build the dynamic matrices on the device
  cudaMalloc((void**) &d_f,ArraySize*sizeof(float));
  cudaMalloc((void**) &d_fnew,ArraySize*sizeof(float));
  //cudaMalloc((void**) &d_matoi,64*sizeof(float));
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] h_f;  delete[] h_fnew;
    cudaFree(d_f);  cudaFree(d_fnew);
    //delete[] h_matoi;  cudaFree(d_matoi);
}
float LatticeBoltzmann::h_rho(int ix,int iy,bool UseNew){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=h_n(ix,iy,i);
    if(UseNew) sum+=h_fnew[n0]; else sum+=h_f[n0];
  }
  return sum;
}
float LatticeBoltzmann::h_Jx(int ix,int iy,bool UseNew,float Fx){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=h_n(ix,iy,i);
    if(UseNew) sum+=h_Vx[i]*h_fnew[n0]; else sum+=h_Vx[i]*h_f[n0];
  }
  return sum+0.5*Fx;
}
float LatticeBoltzmann::h_Jy(int ix,int iy,bool UseNew,float Fy){
  float sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=h_n(ix,iy,i);
    if(UseNew) sum+=h_Vy[i]*h_fnew[n0]; else sum+=h_Vy[i]*h_f[n0];
  }
  return sum+0.5*Fy;
}
float  LatticeBoltzmann::h_Fi(float Ux0,float Uy0,float Fx,float Fy,int i){
  float UdotVi=Ux0*h_Vx[i]+Uy0*h_Vy[i];
  float FdotVi=Fx*h_Vx[i]+Fy*h_Vy[i], UdotF=Ux0*Fx+Uy0*Fy;
  return h_tau[3]*h_w[i]*(FdotVi-UdotF+3*UdotVi*FdotVi);
}
float  LatticeBoltzmann::h_feq(float rho0,float Ux0,float Uy0,int i){
  float UdotVi=Ux0*h_Vx[i]+Uy0*h_Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0*h_w[i]*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}
void LatticeBoltzmann::Start(float rho0,float Ux0,float Uy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
        n0=h_n(ix,iy,i);
        h_f[n0]=h_feq(rho0,Ux0,Uy0,i);
      }
  //Send to the Device
  cudaMemcpy(d_f,h_f,ArraySize*sizeof(float),cudaMemcpyHostToDevice);
}
void LatticeBoltzmann::Collision(float gx,float gy){
    //Do everything on the Device
    dim3 ThreadsPerBlock(N,1,1);
    dim3 BlocksPerGrid(M,1,1);
    d_Collision<<<BlocksPerGrid,ThreadsPerBlock>>>(gx,gy,d_f,d_fnew);
}
void LatticeBoltzmann::ImposeFields(void){
  int iy; float rho_bot,rho_top;
  //Rho for lower wall
  iy=0;
  rho_bot=h_rho(0,iy,false);
  //Rho for upper wall
  iy=Ly-1;
  rho_top=h_rho(0,iy,false);
  //Tarea para el device
  dim3 ThreadsPerBlock(1,1,1); //One thread is simpler
  dim3 BlocksPerGrid(1,1,1);
  d_ImposeFields<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew,rho_bot,rho_top);
  }
void LatticeBoltzmann::Advection(void){
  //Do everything on the Device
  dim3 ThreadsPerBlock(N,1,1);
  dim3 BlocksPerGrid(M,1,1);
  d_Advection<<<BlocksPerGrid,ThreadsPerBlock>>>(d_f,d_fnew);
}
void LatticeBoltzmann::Print(const char * NameFile,float gx,float gy){
  ofstream MyFile(NameFile); float rho0,Ux0,Uy0; int ix,iy; float Fx,Fy;
  //Bring back the data from Device to Host
  //cudaMemcpy(h_f,d_f,ArraySize*sizeof(float),cudaMemcpyDeviceToHost);
  cudaMemcpy(h_fnew,d_fnew,ArraySize*sizeof(float),cudaMemcpyDeviceToHost);
  //Print
  ix=0;
    for(iy=0;iy<Ly;iy++){
      rho0=h_rho(ix,iy,true); Fx=gx*rho0; Fy=gy*rho0;
      Ux0=h_Jx(ix,iy,true,Fx)/rho0; Uy0=h_Jy(ix,iy,true,Fy)/rho0;
      //MyFile<<rho0<<endl;

      MyFile<<iy<<" "<<Ux0<<endl;
  }
  MyFile.close();
}
/*void LatticeBoltzmann::Print(const char * NameFile,float gx,float gy){
  ofstream MyFile(NameFile);
  //Bring back the data from Device to Host
  cudaMemcpy(h_matoi,d_matoi,ArraySize*sizeof(float),cudaMemcpyDeviceToHost);
  //Print
  int ci;
  for(ci=0;ci<64;ci++){
    MyFile<<h_matoi[ci]<<endl;
  }
  MyFile.close();
}*/
//------------ Global main function -------------
int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=100000;
  float rho0=1.0,g=0.01;

  //Start
  Aire.Start(rho0,0,0);
  //Run
  for(t=0;t<tmax;t++){
    Aire.Collision(g,0);
    Aire.ImposeFields();
    Aire.Advection();
  }
  //print
  Aire.Print("Poiseuille_CUDA.dat",g,0);

  return 0;
}
