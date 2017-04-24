#include <iostream>
#include <fstream>
#include <climits>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <algorithm>


//to determine the reaction type
void SetRT(char* rt,int L,int I1A, int I1D, int sTSS, int sTES){
  int i,j,r;

  //A 0 L-1: sX elongation i>i+1
  for (r=0; r<(L-1);r++){
    if(r<sTES){
      rt[r]=2; 
    }//if
    else{
      rt[r]=0; 
    }//else: beyond sense pA site, so no more elongation allowed
  }//for r 

  for (i=0; i<L;i++){
    r=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    if(i>=(sTSS+I1A)){
      rt[r]=1;
    }//if
    else{
      rt[r]=0;
    }//else

    for (j=0; j<L;j++){
      r=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      if(i<(sTSS+I1A) && i>=I1D){
	rt[r]=1;
      }//if
      else{
	rt[r]=0;
      }//else
    
      r=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
      if(i<(sTSS+I1A) && i>=I1D){
	rt[r]=1;
      }//if
      else{
	rt[r]=0;
      }//else
    }//for j
  }//for i

  r=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
  rt[r]=3;
  r=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  rt[r]=4;
  r=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  rt[r]=1;
  r=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  rt[r]=0;//1;!!!!!!!!!!!!!This reaction is not allowed anymore!
  r=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 
  rt[r]=1;
  r=2*L*L+2*L+4;                //J 7 1: sFLC > void
  rt[r]=1;
  r=2*L*L+2*L+5;                //K 8 1: OFF > ON
  rt[r]=5;
  r=2*L*L+2*L+6;                //L 9 1: ON > OFF
  rt[r]=1;
}//SetRT


//k0 = PolII elongation rate
//k1 = splicing
//k2 = 53 intronic RNA degradation rate
//k3 = 35 intronic RNA degradation rate
//k4 = P2 drop off (after reaching pA site)
//k5 = sP initiation rate
//k6 = sense FLC export rate
//k7 = sFLC degradation rate
//k8 = OFF>ON rate
//k9 = ON>OFF rate
//array with rates
void SetK(double ** k,int L,double k0,double k1,double k2,double k3,double k4,double k5,double k6,double k7,double k8, double k9){
  int r;

  //A 0 L-1: sX elongation i>i+1
  for (r=0; r<(L-1);r++){
    k[r][0]=k0;
  }//for r 
  //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
  for (r=(L-1); r<(2*L-1);r++){
    k[r][0]=k1;
  }//for r
  //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
  for (r=2*L-1; r<(L*L+2*L-1);r++){
    k[r][0]=k2; 
  }//for r
  //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
  for (r=(L*L+2*L-1); r<(2*L*L+2*L-1);r++){
    k[r][0]=k3; 
  }//for r
  r=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
  k[r][0]=k4;
  r=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  k[r][0]=k5;
  r=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  k[r][0]=k1;
  r=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  k[r][0]=k6;
  r=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 
  k[r][0]=k6;
  r=2*L*L+2*L+4;                //J 7 1: sFLC > void
  k[r][0]=k7;
  r=2*L*L+2*L+5;                //K 8 1: OFF > ON
  k[r][0]=k8;
  r=2*L*L+2*L+6;                //L 9 1: ON > OFF
  k[r][0]=k9;
 
}//SetK


void SetR_membernumber(int* m, int L){
  int i,j,r;

  for(i=0;i<L-1;i++){ 
    r=i;                    //A 0 L-1: sX elongation i>i+1
    m[r]=4;    
  }//for i <L-1

  for(i=0; i<L;i++){
    for(j=0; j<L;j++){
      r=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      m[r]=1;
      r=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
      m[r]=1;
    }//for j
    r=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    m[r]=1;
  }//for i < L
  r=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
  m[r]=2;
  r=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  m[r]=2;
  r=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  m[r]=1;
  r=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  m[r]=1;
  r=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 
  m[r]=1;
  r=2*L*L+2*L+4;                //J 7 1: sFLC > void
  m[r]=1;
  r=2*L*L+2*L+5;                //K 8 1: OFF > ON
  m[r]=1;
  r=2*L*L+2*L+6;                //L 9 1: ON > OFF
  m[r]=1;
}//SetR_membernumber


void SetR(int*** R,int L,int* sP,int* sS1,int** IN,int* sFLC,int* unsFLC,int* s1FLC,int* STATE,int sTSS, int sTES){
  int r,i,j;

 for(i=0;i<L-1;i++){ 
    r=i;                    //A 0 L-1: sX elongation i>i+1
    R[r][0]=sP+i;
    R[r][1]=sS1+i;
    R[r][2]=sP+i+1;
    R[r][3]=sS1+i+1; 
  }//for i <L-1

  for(i=0; i<L;i++){
    for(j=0; j<L;j++){
      r=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      R[r][0]=IN[i]+j;
      r=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
      R[r][0]=IN[i]+j;
    }//for j
    r=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    R[r][0]=sP+i;
  }//for i < L
  r=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
  R[r][0]=sP+sTES;
  R[r][1]=sS1+sTES;
  r=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  R[r][0]=sP+sTSS; 
  R[r][1]=STATE;
  r=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  R[r][0]=unsFLC;
  r=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  R[r][0]=unsFLC;
  r=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 
  R[r][0]=s1FLC;
  r=2*L*L+2*L+4;                //J 7 1: sFLC > void
  R[r][0]=sFLC;
  r=2*L*L+2*L+5;                //K 8 1: OFF > ON
  R[r][0]=STATE;
  r=2*L*L+2*L+6;                //L 9 1: ON > OFF
  R[r][0]=STATE;
}//SetR


void SetSizecon (int *size,int L,int sTSS,int sTES){
  int i,j,r;//r reaction index
  
  for(i=0;i<L-1;i++){ 
    r=i;                    //A 0 L-1: sX elongation i>i+1
    size[r]=5;
    if(i==sTSS){
      size[r]=6;     
    }//if 
    else if(i==(sTES-1)){
      size[r]=6;
    }//else if 
  }//for i <L-1

  for(i=0; i<L;i++){
    for(j=0; j<L;j++){
      r=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      size[r]=4;
      r=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
      size[r]=4;
    }//for j
    r=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    size[r]=3;
  }//for i < L

  r=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
  size[r]=7;
  r=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  size[r]=2;
  r=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  size[r]=5;
  r=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  size[r]=5;
  r=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 
  size[r]=2;
  r=2*L*L+2*L+4;                //J 7 1: sFLC > void
  size[r]=1;
  r=2*L*L+2*L+5;                //K 8 1: OFF > ON
  size[r]=3;
  r=2*L*L+2*L+6;                //L 9 1: ON > OFF
  size[r]=3;
}//SetSizecon

//from reaction index r to site i
//i=r;                                 //A 0 L-1: sX elongation i>i+1 
//i=r-(L-1);                           //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
//i=(r-2*L+1)%L; j=(r-2*L+1)/L;        //C 2 L*L: IN[i,j] > IN[i+1,j] (void if i+1=j)
//i=(r-L*L-2*L+1)%L; j=(r-L*L-2*L+1)/L;//D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
//na                                   //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
//na                                   //F 5 1: void > sP[sTSS]
//na                                   //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
//na                                   //H 6 1: unsFLC > sFLC +IN[I1D,I1A]
//na                                   //I 6 1: s1FLC > sFLC 
//na                                   //J 7 1: sFLC > void
//na                                   //K 8 1: OFF > ON
//na                                   //L 9 1: ON > OFF

//Vice Versa from site i to reaction index r, site i:
//r=i;                    //A 0 L-1: sX elongation i>i+1 
//r=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
//r=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
//r=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
//r=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
//r=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
//r=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
//r=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
//r=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 
//r=2*L*L+2*L+4;                //J 7 1: sFLC > void
//r=2*L*L+2*L+5;                //K 8 1: OFF > ON
//r=2*L*L+2*L+6;                //L 9 1: ON > OFF


void SetCon (int ** con,int L,int sTSS,int sTES,int I1D,int I1A){
  int i,j,r;//i coord of site, r reaction index

  for(i=0;i<L-1;i++){ 
    r=i;                    //A 0 L-1: sX elongation i>i+1
    con[r][0]=i-1;                    //A 0 L-1: sX elongation i>i+1
    con[r][1]=i;                      //A 0 L-1: sX elongation i>i+1
    con[r][2]=i+1;                    //A 0 L-1: sX elongation i>i+1
    con[r][3]=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    con[r][4]=i+1+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    if(i==sTSS){
      con[r][5]=2*L*L+2*L;                  //F 5 1: void > sP[sTSS]
    }//if 
    else if(i==(sTES-1)){
      con[r][5]=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
    }//else if 
  }//for i <L-1

  for(i=0; i<L;i++){
    for(j=0; j<L;j++){
      r=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      con[r][0]=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      con[r][1]=i+1+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      con[r][2]=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
      con[r][3]=i+1+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)

      r=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
      con[r][0]=i+j*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      con[r][1]=i+(j-1)*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
      con[r][2]=i+j*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
      con[r][3]=i+(j-1)*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
    }//for j

    r=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    con[r][0]=i+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
    con[r][1]=I1D+(sTSS+I1A-1)*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
    con[r][2]=I1D+(sTSS+I1A-1)*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)

  }//for i < L

  r=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
  con[r][0]=sTES-1;                    //A 0 L-1: sX elongation i>i+1
  con[r][1]=sTES;                    //A 0 L-1: sX elongation i>i+1
  con[r][2]=sTES+(L-1);              //B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
  con[r][3]=2*L*L+2*L-1;          //E 4 1: sX[sTES] > unsFLC (or s1sFLC)
  con[r][4]=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  con[r][5]=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  con[r][6]=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 

  r=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  con[r][0]=sTSS;                    //A 0 L-1: sX elongation i>i+1
  con[r][1]=2*L*L+2*L;               //F 5 1: void (+STATE) > sP[sTSS]

  r=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  con[r][0]=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  con[r][1]=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  con[r][2]=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC
  con[r][3]=I1D+(sTSS+I1A-1)*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
  con[r][4]=I1D+(sTSS+I1A-1)*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)

  r=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  con[r][0]=2*L*L+2*L+1;                //G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  con[r][1]=2*L*L+2*L+2;                //H 6 1: unsFLC > sFLC + IN[I1D,I1A]
  con[r][2]=2*L*L+2*L+4;                //J 7 1: sFLC > void
  con[r][3]=I1D+(sTSS+I1A-1)*L+2*L-1;          //C 2 L*L IN[i,j] > IN[i+1,j] (void if i+1=j)
  con[r][4]=I1D+(sTSS+I1A-1)*L+L*L+2*L-1;      //D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)

  r=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC 
  con[r][0]=2*L*L+2*L+3;                //I 6 1: s1FLC > sFLC
  con[r][1]=2*L*L+2*L+4;                //J 7 1: sFLC > void
  
  r=2*L*L+2*L+4;                //J 7 1: sFLC > void
  con[r][0]=2*L*L+2*L+4;                //J 7 1: sFLC > void
 
  r=2*L*L+2*L+5;                //K 8 1: OFF > ON
  con[r][0]=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  con[r][1]=2*L*L+2*L+5;                //K 8 1: OFF > ON
  con[r][2]=2*L*L+2*L+6;                //L 9 1: ON > OFF

  r=2*L*L+2*L+6;                //L 9 1: ON > OFF
  con[r][0]=2*L*L+2*L;                  //F 5 1: void (+STATE) > sP[sTSS]
  con[r][1]=2*L*L+2*L+5;                //K 8 1: OFF > ON
  con[r][2]=2*L*L+2*L+6;                //L 9 1: ON > OFF

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
