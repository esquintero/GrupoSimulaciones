#include <iostream>
#include <string>
//#include <ifstream>
using namespace std;
int main(void){
  //string filename;
  //cout<<"Please enter the name of the .dat file: ";
  //cin>> filename;
  //cout<<"set size square "<<endl;
  cout<<"set multiplot"<<endl;
  cout<<"set grid"<<endl;
  //cout<<"plot \""<<filename<<"\" using 1:2 w l"<<endl;
  cout<<"plot \"S-J_original_axis.dat\" using 3:4 title \"Jupiter\" w l,";
  cout<<" \"S-J_original_axis.dat\" using 1:2 title \"Sol\" w l,";
  cout<<" \"S-J_original_axis.dat\" using 5:6 title \"Troyano\" w l"<<endl;
  cout<<"set origin .25,.25"<<endl;
  cout<<"set size square .5,.5"<<endl;
  cout<<"set label \"Sol\""<<endl;
  cout<<"clear"<<endl;
  cout<<"unset key"<<endl;
  cout<<"unset object"<<endl;
  cout<<"plot \"S-J_original_axis.dat\" using 1:2 title \"Sol\" w l linecolor \"dark-green\""<<endl;
  cout<<"unset multiplot"<<endl;
  cout<<"pause 10"<<endl;

  return 0;
}
