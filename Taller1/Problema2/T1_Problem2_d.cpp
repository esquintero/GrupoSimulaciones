#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;
const double ErrMax=1e-7;

double f(double alpha, double r, double l, double t){
  return cos(alpha*t-l*r*sin(t));
}

double IntegralPorSimpson(double alpha, double r, double l, double a, double b, int n){
    double t,h,sum; int i;
    n=2*n; 
    h=(b-a)/n;
    sum=0;
    for(i=0;i<=n;i++){
        t=a+i*h;
        if(i==0 || i==n)
            sum=sum+f(alpha,r,l,t);
        
        else if(i%2==0)
            sum=sum+2*f(alpha,r,l,t);
        
        else
            sum=sum+4*f(alpha,r,l,t);
        
    }
    return sum*h/3;
}

double Bessel(double alpha,double r, double l){
    double a=0,b=M_PI; int n=50;
    return 1.0/M_PI*IntegralPorSimpson(alpha,r,l,a,b,n);
}

double ZerosBi(double alpha,double r, double l, double a, double b){
    double x;
    double m, fa, fm;

    fa=Bessel(alpha,r,a);
    while(b-a >= ErrMax){
        m = (b+a)/2; fm=Bessel(alpha,r,m);
        if(fa*fm>0)
        {a=m; fa=fm;}
        else
        b=m;
        }
        return (a+b)/2;
}

int main(){
    double alpha=0;
    double r=1; 
    double l;
    cout<<"El primer cero está en "<<ZerosBi(alpha,r,l,0.01,4)<<endl;
    cout<<"El segundo cero está en "<<ZerosBi(alpha,r,l,4,8)<<endl;
    cout<<"El tercer cero está en "<<ZerosBi(alpha,r,l,8,12)<<endl;
    cout<<"El cuarto cero está en "<<ZerosBi(alpha,r,l,12,15)<<endl;
}
