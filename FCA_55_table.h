#include <iostream>
#include <fstream>
#include <climits>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>


//to determine the reaction type
void SetRT(char* rt){
  int r;
  r=0;                   //A 0 1: void > FLC
  rt[r]=0;
  r=1;                   //B 0 1: FLC > void
  rt[r]=1;
  r=2;                   //C 0 1: void > void
  rt[r]=3;
  r=3;                //D 8 1: OFF > ON
  rt[r]=2;
  r=4;                //E 9 1: ON > OFF
  rt[r]=1;
}//SetRT


//k0 = Cellular FLC production rate
//k1 = FLC degradation rate
//k2 OFF > ON rate
//k3 ON > OFF rate
//array with rates
void SetK(double * k,double k0,double k1, double k2, double k3){
  int r;
  r=0;                   //A 0 1: void > FLC
  k[r]=k0;
  r=1;                   //B 1 1: FLC > void
  k[r]=k1;
  r=2;                   //C 1 1: void > void
  k[r]=k3;
  r=3;                //D 2 1: OFF > ON
  k[r]=k2;
  r=4;                //E 3 1: ON > OFF
  k[r]=k3;
 
}//SetK


void SetR_membernumber(int* m){
  int r;
  r=0;                   //A 0 1: void > FLC
  m[r]=1;
  r=1;                   //B 1 1: FLC > void
  m[r]=1;
  r=2;                   //C 1 1: void > void
  m[r]=1;
  r=3;                   //D 2 1: OFF > ON
  m[r]=1;
  r=4;                   //E 3 1: ON > OFF
  m[r]=1;
}//SetR_membernumber


void SetR(int*** R,int* sFLC, int* STATE){
  int r;
  r=0;                   //A 0 1: void > FLC
  R[r][0]=STATE;
  r=1;                   //B 1 1: FLC > void
  R[r][0]=sFLC;
  r=2;                   //C 1 1: void > void
  R[r][0]=sFLC;
  r=3;                   //D 2 1: OFF > ON
  R[r][0]=STATE;
  r=4;                   //E 3 1: ON > OFF
  R[r][0]=STATE;
}//SetR


void SetSizecon (int *size){
  int r;//r reaction index
  r=0;                   //A 0 1: void > FLC
  size[r]=2;
  r=1;                   //B 0 1: FLC > void
  size[r]=1;
  r=2;                   //C 1 1: void > void
  size[r]=1;
  r=3;                   //D 2 1: OFF > ON
  size[r]=3;
  r=4;                   //E 3 1: ON > OFF
  size[r]=3;
}//SetSizecon


//r = reaction index
//r=0;               //A 0 1: void > FLC
//r=1;               //B 0 1: FLC > void
//r=2;                   //C 1 1: void > void
//r=3;                   //D 2 1: OFF > ON
//r=4;                   //E 3 1: ON > OFF
void SetCon (int ** con){
  int r;//r reaction index
    r=0;                   //A 0 1: void > FLC
    con[r][0]=0;           //B 0 1: FLC > void
    con[r][0]=1;           //B 0 1: FLC > void
    r=1;                   //B 0 1: FLC > void
    con[r][0]=1;           //B 0 1: FLC > void
    r=2;                   //C 1 1: void > void
    con[r][0]=2;           //C 1 1: void > void
    r=3;                   //D 2 1: OFF > ON
    con[r][0]=0;           //A 0 1: void > FLC
    con[r][1]=3;           //D 2 1: OFF > ON
    con[r][2]=4;           //E 3 1: ON > OFF
    r=4;                   //E 3 1: ON > OFF
    con[r][0]=0;           //A 0 1: void > FLC
    con[r][1]=3;           //D 2 1: OFF > ON
    con[r][2]=4;           //E 3 1: ON > OFF
}//SetCon


//store this part in here, might need it later on

//   //  Temp test
//   cout << endl;
//   cout << endl << "RT" << endl;
//   for(i=0;i<NR;i++){
//     cout << "in"<< i << " ";
//     if(reactiontype[i]==1){
//       cout  << "1 ";
//     }
//     else if(reactiontype[i]==2){
//       cout  << "2 ";
//     }
//     else if(reactiontype[i]==3){
//       cout  << "3 ";
//     }
//     else if(reactiontype[i]==4){
//       cout  << "4 ";
//     }
//     else if(reactiontype[i]==5){
//       cout  << "5 ";
//     }
//     else if(reactiontype[i]==6){
//       cout  << "6 ";
//     }
//     else if(reactiontype[i]==7){
//       cout  << "7 ";
//     }
//     else if(reactiontype[i]==0){
//       cout  << "0 ";
//     }  
//   cout << endl;
//   }//for i
//   cout << endl << "k" << endl;
//   for(i=0;i<NR;i++){
//     cout << i << " " << k[i][0] << " " << k[i][1]  << " " << k[i][2]  << endl;;
//   }//for i

//   cout << endl << "r[0]" << endl;  
//   for(i=0;i<NR;i++){
//     cout << i << " "<< *r[i][0] << " ";
    
//   }
//   cout << endl;
//   cout << endl << "r[1]" << endl;
//   for(i=0;i<NR;i++){
//     if(r_membernumber[i]>1){
//       cout << i << " "<< *r[i][1] << " ";
//     }
//   }
//   cout << endl;
//   cout << endl << "connect " << endl;
  
//   for(i=0;i<NR;i++){
//     cout << i << " ";
//     for (j=0;j<sizecon[i];j++){
//       cout << connect[i][j] << " " ;
//     }
//     cout << endl;
//   }
//   cout << endl << "p" << endl;
//   for(i=0;i<NR;i++){
//     cout << i << " "<< p[i] << " ";
//   }

//   cout << endl;
//  //TEMp
