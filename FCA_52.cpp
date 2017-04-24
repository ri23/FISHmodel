//Robert Ietswaart
//20160924
//FCA_52.cpp

#include <iostream>
#include <fstream>
#include <math.h>
#include <ctime>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <cstring>
#include "MersenneTwister.h"
#include "FCA_52_table.h"

using namespace std;

const double Time = 10*24*60*60;//Time in s
const int cap =1000;//artificial cap that HMT  and FLD levels shouldn't exceed
const int L=209,I1A=127,I1D=14,sTSS=4,asTSS=208,sTES=204;// 30 bp(Pol II footprint),R2= region to switch to S2 (see review)

MTRand mtrand1;       //    = properly seeding
//MTRand mtrand1(1973);
//MTRand mtrand1(1000);   // for test, so that seed=1973 is the same for all testruns, to detect bug more easily!

inline char* itoa(long int i){
  static char a[256]={0};
  sprintf(a,"%ld",i);
  return a;
}//convert from int to string (char array)


int GetGroupNumber (double* p_g, double p_s, int N_g, bool& faal){
  double help=0.0;
  double r2= mtrand1.rand(p_s);//uniform in [0,p_s]
  int i, groupnumber=-1;
  for(i=1;i<N_g;i++){
    help=help+p_g[i];
    if(r2<=help){
      groupnumber = i;
      i=N_g-1;//we are done so set i to end value
    }//if
  }//for i
  if( (groupnumber <= 0) || (groupnumber >= N_g)){//mag niet, er wordt gesampled uit 1..N_g
    cout << "faalaap" << endl;
    faal = 1;
  }//if
  return groupnumber;
}//GetGroupNumber


//membernumber = aantal members in group groupnumber (>=1, <NR)
//groupnumber is index of group (>=1)
int RejectionSample ( int** index, double * p,  int membernumber, int groupnumber, double pmin,  bool & faal, double t){
  int reaction;
  double r3;
  int n=0,r4;

  do {
    r3=pmin*pow(2.0,double(groupnumber));
    r3= mtrand1.rand( r3 );
    r4=membernumber-1;
    r4=mtrand1.randInt( r4 ); //should be int in 0..(membernumber-1)
    reaction= index[groupnumber][r4];
 
    if (r4 < membernumber){
      reaction= index[groupnumber][r4];
    }//if
    else{
      cout << "fail,groupnumber=" << groupnumber<< " membernumber=" << membernumber << " x=" << r4 << " reaction=" << reaction << " p[reaction]=" << p[reaction] << endl;
      reaction= index[groupnumber][(r4-1)]; //dit is natuurlijk een faaloplossing, met MT moet dit gewoon goed gaan altijd!
    }//else

    //temp
    n++;
    if (n>=100){
      cout << "inf sampling loop, groupnumber=" << groupnumber<< " membernumber=" << membernumber << " x= " << r4 << " reaction=" << reaction << " p[reaction]=" << p[reaction] << endl;
      faal = 1;
      break;
    }//if checktime
    //temp

  }while(r3 > p[reaction]);
  return reaction;
}//Sample

double SetPmin (double k0,double k1,double k2,double k3,double k4,double k5,double k6,double k7,double k8,double k9, int& Ng){
  int arraysize=9;//arraysize is number of rates involved
  double kmin, k_sort[arraysize];
  k_sort[0]=k0;
  k_sort[1]=k1;
  k_sort[2]=k2;
  k_sort[3]=k3;
  k_sort[4]=k4;
  k_sort[5]=k5;
  k_sort[6]=k6;
  k_sort[7]=k7;
  k_sort[8]=k8;
  // k_sort[9]=k9;

  sort(k_sort,k_sort+arraysize);
  kmin= k_sort[0];

  k_sort[0]=k0;
  k_sort[1]=k1*cap;//unsFLC
  k_sort[2]=k2*cap;
  k_sort[3]=k3*cap;
  k_sort[4]=k4;
  k_sort[5]=k5;
  k_sort[6]=k6*cap;
  k_sort[7]=k7*cap;
  k_sort[8]=k8;
  //k_sort[9]=k9;

  sort(k_sort,k_sort+arraysize);
  cout << "pmax =  " << k_sort[arraysize-1] << endl;
  //Ng is smallest integer such that pmax < pmin*2^Ng

  frexp((k_sort[arraysize-1]/kmin),&Ng);
  Ng++;//1 bigger because you also have group 0
  return kmin;
}//SetPmin


//known site has to be deleted from group
void DeleteReactionToGroup (int** group2reaction,int site, int* membernumber, int** reaction2group){
  int group, index;
  group = reaction2group[site][0];
  index = reaction2group[site][1]; 
  group2reaction[group][index]= group2reaction[group][(membernumber[group]-1)];//last site one goes into place of to be deleted one
  reaction2group[group2reaction[group][(membernumber[group]-1)] ][0]= group;// fix reaction2group of last site, not deleted one!
  reaction2group[group2reaction[group][(membernumber[group]-1)] ][1]= index;
  group2reaction[group][(membernumber[group]-1)]=-1;//delete last element of group2reaction
  membernumber[group]=membernumber[group]-1;//now the group has one member less
  reaction2group[site][0]= -1;
  reaction2group[site][1]= -1;
}//DeleteReactionFromIndex

//NB if propensity is zero. it is in group=0 
//always first delete from old group/index, then add to new one 
//site is a known reaction index, groupnumber calculated via frexp
void AddReactionToGroup(int** group2reaction, int group, int site, int* membernumber, int** reaction2group){
  group2reaction[group][membernumber[group]]=site;//added at end of group2reaction
  reaction2group[site][0]=group;//fix reaction2group
  reaction2group[site][1]=membernumber[group];
  membernumber[group]=membernumber[group]+1;
}//AddReactionToGroup

double ElonXP(double* k, int** r){
  double pnew;
  if((*r[2]+*r[3])==0){
    pnew=(*r[0]+*r[1])*k[0];
  }//if new site not occupied by P2
  else{
    pnew=0;
  }//if new site occupied, can't elongate
  //r[0]=sP at site;
  //r[1]=sS1 at site;
  //r[2]=sP at neighboring site
  //r[3]=sS1 at neighboring site
  return pnew;
}//ElonXP
 
double LoadP(double* k, int** r){
  double pnew;
  if((*r[0]==0) && (*r[1]==1)){
    pnew=k[0];
  }//if new site not occupied by P2
  else{
    pnew=0;
  }//if new site occupied, can't fire
  //r[0]=sP at TSS site
  //r[1]= STATE
  return pnew;
}//LoadElonP


//site is reaction index
double CalculateP (double** k, int*** r, char RT, int site){
  double pnew;
  if (RT==0){
    pnew=0;
  }//virtual reaction
  else if (RT==1){
    pnew = k[site][0]*(*r[site][0]);//first order reaction   
  }//first order reaction (B CoT splicing, C D IN[i,j] degradation, G H I J RNA metabolism, L switch to OFF STATE
  else if (RT==2){
    pnew = ElonXP(k[site],r[site]);
  }// if A sX elongation
  else if (RT==3){
    pnew= k[site][0]*(*r[site][0]+*r[site][1]);//sP+sS1 
  }//if E  Pol II drop off
  else if (RT==4){
    pnew= LoadP(k[site],r[site]);
  }//if F sP firing
  else if (RT==5){    
    pnew = k[site][0]*(1-*r[site][0]);
  }// if K switch to ON STATE
  else{
    cout << "faal, dit reactietype is niet toegestaan " << RT << " " << k[site] << " " << k[site] << " site="<< site << endl;
    pnew=0;//temp faal oplossing, deze shit mag eenvoudig niet voorkomen
  }// else
  // cout << pnew << endl;
  return pnew;
}//double CalculateP


void SetP (double * p, double** k, int*** r, char* RT,  int * g_membernumber,int ** group2reaction, int NR ,double pmin, int** reaction2group){
  int i, g;
  for(i=0;i< NR; i++){// initial propensity reactions
    // cout << "SetP0 i=" << i << endl;
    p[i]=CalculateP(k,r,RT[i],i);
    frexp((p[i]/pmin),&g);
    // cout << "SetP1 p[i]=" << p[i] << " p[i]/pmin" << p[i]/pmin <<  endl;
    AddReactionToGroup(group2reaction,g,i,g_membernumber, reaction2group);
    //  cout << "SetP2 g=" << g << endl << endl;
  }//for i
}//SetP


void SetPgroup (double * p, double * p_group, int N_group, int * membernumber, int ** index){
  int i,j;
  for(i=0;i<N_group;i++){
    p_group[i]=0;
    for(j=0;j<membernumber[i];j++){
      p_group[i]= p_group[i]+p[index[i][j]];
    }//for j
  }//for
}//SetPgroup


double SetPsum ( double* p_g, int N_g){
  int i;
  double p_s=0; 
  for(i=1;i<N_g;i++){//only nonzero propensity groups (1..N_g) are added
    p_s=p_s+p_g[i];
  }//for
  return p_s;
}//setPsum


void UpdatePgroup(double pold, double pnew, double* pg, int group, double pmin){
  pg[group]=pg[group]+pnew-pold;
  if (pg[group]<((pmin*8)/10)){//Shabby solution for problem of floating point values
    pg[group]=0;
  }//if
}//UpdatePgroup


//sum over g=1..Ng-1, the groups with nonzero propensity
double UpdatePsum(double* pg, int Ng){
  int i;
  double help=0;
  for(i=1;i<Ng; i++){
    help = help + pg[i];
  }//for i 
  return help;
}//UpdatePsum


void UpdateP2 (double* p, double pnew, int** group2reaction, int* g_membernumber, int** reaction2group, double* pg, double pmin, int site){
  int g;//groupnumber
  if (pnew != p[site]){
    frexp((pnew/pmin),&g);//calculate new groupnumber
    //cout << "Up2 0 g="<< g << endl;
    if (g!= (reaction2group[site][0])){
      UpdatePgroup(p[site],0,pg, reaction2group[site][0],pmin);//old group
      UpdatePgroup(0,pnew, pg,g,pmin);//new group
      DeleteReactionToGroup (group2reaction,site ,g_membernumber,reaction2group);
      AddReactionToGroup(group2reaction, g,site,g_membernumber, reaction2group);
      p[site]=pnew;  
      //cout << "Up2 1 p=" << p[site] << " site=" << site << endl;
    }//if different group
    else{
      UpdatePgroup(p[site],pnew,pg, reaction2group[site][0],pmin);//same group
      p[site]=pnew;  
      //cout << "Up2 2 p=" << p[site] <<" site=" << site << endl;
    }//pnew still in same group
  }//if pnew!= p[site], anders hoef je geen actie te ondernemen
}//UpdateP2


//m is current reaction
void PerformReaction(int m, int*** r, int** connect, int* sizecon,int* sP,int* sS1,int** IN,int& sFLC,int& unsFLC,int& s1FLC,int& STATE,double t,bool& fout){

  int i,j;//positions on grid

  if( m < (L-1) ){
    i=m;
    if (sP[i]>0){
      sP[i]--;
      sP[i+1]++;
    }//if sP
    else if (sS1[i]>0){
      sS1[i]--;
      sS1[i+1]++;
    }//if sS1
    //temp
    if (sP[i]<0 || sP[i+1]>1 || sS1[i]<0 || sS1[i+1]>1){
      cout << endl << "cA " << "sP[i]= " << sP[i] << "sP[i+1]= " << sP[i+1] <<"sS1[i]= " << sS1[i] << "sS1[i+1]= " << sS1[i+1] << " m=" << m << " i="<< i <<  endl; 
      fout=1;
    }
    //temp
  }//A 0 L-1: sX elongation i>i+1
  else if(m < (2*L-1)){
    i=m-(L-1);
    sP[i]--;
    sS1[i]++;
    IN[I1D][sTSS+I1A-1]++;
    //temp
    if (sP[i]<0 || sS1[i]>1 || IN[I1D][sTSS+I1A-1]>cap ){
      cout << endl << "cB " << "sP[i]= " << sP[i] << " sS1[i]= " << sS1[i] << "IN[I1D][I1A]= " << IN[I1D][sTSS+I1A-1] << " m=" << m << " i=" << i << endl; 
      fout=1;
    }
    //temp 
  }//B 1 L: sP[i] > sS1[i] + IN[I1D,I1A]
  else if (m <(L*L+2*L-1)){
    i=(m-2*L+1)%L; 
    j=(m-2*L+1)/L;       
    IN[i][j]--;
    IN[i+1][j]++;
    if(i+1==j){
      IN[i+1][j]--;
    }//if end of transcript, then reaction: IN > void
    //temp
    if (IN[i][j]<0 || IN[i+1][j]>cap){
      cout << endl << "cC =" << "IN[i][j]= " << IN[i][j] <<  "IN[i+1][j]= " << IN[i+1][j] << " m=" << m << " i=" << i << " j=" << j << endl; 
      fout=1;
    }
    //temp
  }//C 2 L*L: IN[i,j] > IN[i+1,j] (void if i+1=j)
  else if (m <(2*L*L+2*L-1)){
    i=(m-L*L-2*L+1)%L; 
    j=(m-L*L-2*L+1)/L;      
    IN[i][j]--;
    IN[i][j-1]++;
    if(i==j-1){
      IN[i][j-1]--;
    }//if end of transcript, then reaction: IN > void
    //temp
    if (IN[i][j]<0 || IN[i][j-1]>cap){
      cout << endl << "cD =" << "IN[i][j]= " << IN[i][j] <<  "IN[i][j-1]= " << IN[i][j-1] << " m=" << m << " i=" << i << " j=" << j << endl; 
      fout=1;
    }
    //temp
  }//D 3 L*L: IN[i,j] > IN[i,j-1] (void if i=j-1)
  else if (m ==(2*L*L+2*L-1) ){          
    if (sP[sTES]>0){
      sP[sTES]--;
      unsFLC++; 
    }//if sP
    else if (sS1[sTES]>0){
      sS1[sTES]--;
      s1FLC++;
    }//if sS1
    //temp
    if (sP[sTES]<0 || unsFLC>cap || sS1[sTES]<0 || s1FLC>cap){
      cout << endl << "cE " << "sP[sTES]= " << sP[sTES] << " sS1[sTES]= " << sS1[sTES] << " unsFLC= " << unsFLC << " s1FLC= " << s1FLC << " m=" << m <<  endl; 
      fout=1;
    }//if
    //temp
  }//E 4 1: sX[sTES] > unsFLC (or s1sFLC) 
  else if (m == (2*L*L+2*L) ){
    sP[sTSS]++;
    //temp
    if (sP[sTSS]>1){
      cout << endl << "cF =" << " sP[TSS]= " << sP[sTSS] << " m=" << m << endl; 
      fout=1;
    }
    //temp
  }//F 5 1: void (+STATE) > sP[sTSS]
  else if (m == (2*L*L+2*L+1) ){
    unsFLC--;
    s1FLC++;
    IN[I1D][sTSS+I1A-1]++;
    //temp
    if (unsFLC<0 || s1FLC>cap || IN[I1D][sTSS+I1A-1]>cap){
      cout << endl << "cG =" << "unsFLC= " << unsFLC << "s1FLC= " << s1FLC << "IN[I1D][I1A]= " << IN[I1D][sTSS+I1A-1] <<  " m=" << m << endl; 
      fout=1;
    }
    //temp
  }//G 1 1: unsFLC  > s1FLC + IN[I1D,I1A]
  else if (m == (2*L*L+2*L+2) ){
    //unsFLC--;
    //sFLC++;
    //IN[I1D][sTSS+I1A-1]++;
    //temp
    //if (unsFLC<0 || IN[I1D][sTSS+I1A-1]>cap || sFLC>cap){
      cout << endl << "cH =" << "unsFLC= " << unsFLC << " IN[I1D][I1A]= " << IN[I1D][sTSS+I1A-1] << " sFLC=" << sFLC << " m=" << m << endl;
      fout=1;
      //}
    //temp
  }//H 6 1: unsFLC > sFLC + IN[I1D,I1A]> NOT ALLOWED ANYMORE, because it speeds up splicing artificially
  else if (m == (2*L*L+2*L+3) ){
    s1FLC--;
    sFLC++;
    //temp
    if (s1FLC<0 || sFLC>cap){
      cout << endl << "cI =" << "s1FLC= " << s1FLC << " sFLC=" << sFLC << " m=" << m << endl;
      fout=1;
    }
    //temp
  }//I 6 1: s1FLC > sFLC
  else if (m == (2*L*L+2*L+4) ){
    sFLC--;
    //temp
    if (sFLC<0){
      cout << endl << "cJ =" << "sFLC= " << sFLC <<  " m=" << m << endl; 
      fout=1;
    }
    //temp
  }//J 7 1: sFLC > void
  else if (m == (2*L*L+2*L+5) ){
    STATE=1;
  }//K 8 1: OFF > ON
  else if (m == (2*L*L+2*L+6) ){
    STATE=0;
  }//L 9 1: ON > OFF
  //temp
  else{
    cout << "faal Perf reacton m=" << m << endl;
  }//else
  //temp
}//PerformReaction

//you only get in here if propensity actually changes
void UpdateP1 (int m, double* p, double** k, int*** r, char* rt, int** connect , int* sizecon,int** group2reaction, int* g_membernumber, int** reaction2group, double* pg, int Ng, double pmin, double& ps){
  int i, site;
  double pnew;
  //  cout << "m= " << m << endl;
  for (i=0; i<(sizecon[m]);i++){
    site=connect[m][i];//site is current reaction dependent on m
    //cout << "site= " << site << endl;
    pnew = CalculateP (k,r,rt[site],site);
    UpdateP2(p,pnew, group2reaction, g_membernumber,  reaction2group,  pg, pmin,site);
  }//for  all dependent reactions
 
  ps=UpdatePsum(pg,Ng);
  if (ps <= 0){
    cout << "3faal psum=" << ps << endl;
  }//if
}//UpdateP

void OutputToFile(int* sP,int* sS1,int** IN,int sFLC,int unsFLC,int s1FLC, double Vol ,string fname){

  ofstream output;
  int i,j,Us=0,Um=0,U5=0,U3=0,Loc=0;
   
  Loc=unsFLC+s1FLC;
  for(i=(sTES-11);i<(sTES+1);i++){//2/3 of mRNA exon signal=11sites upstream of sTES
    Loc=Loc+sP[i]+sS1[i];
  }//for i
  Us=unsFLC;
  for(i=I1D+(sTSS+I1A-I1D)*2/3;i<L;i++){
    Us=Us+sP[i];
  }//for i
  for(i=I1D;i<I1D+(sTSS+I1A-I1D)/3;i++){
    for(j=i+(sTSS+I1A-I1D)*2/3;j<(sTSS+I1A);j++){
      Us=Us+IN[i][j];
    }//for j
  }//for i
  Um=unsFLC;
  for(i=I1D+(sTSS+I1A-I1D)*7/12;i<L;i++){
    Um=Um+sP[i];
  }//for i
  for(i=I1D;i<I1D+(sTSS+I1A-I1D)/4;i++){
    for(j=I1D+(sTSS+I1A-I1D)*7/12;j<(sTSS+I1A);j++){
      Um=Um+IN[i][j];
    }//for j
  }//for i
  for(i=I1D+(sTSS+I1A-I1D)/4;i<I1D+(sTSS+I1A-I1D)*5/12;i++){
    for(j=i+(sTSS+I1A-I1D)/3;j<(sTSS+I1A);j++){
      Um=Um+IN[i][j];
    }//for j
  }//for i
  U5=unsFLC;
  for(i=I1D+(sTSS+I1A-I1D)/3;i<L;i++){
    U5=U5+sP[i];
  }//for i
  for(i=I1D;i<I1D+(sTSS+I1A-I1D)/6;i++){
    for(j=i+(sTSS+I1A-I1D)/3;j<I1D+(sTSS+I1A-I1D)/2;j++){
      U5=U5+IN[i][j];
    }//for j
    for(j=I1D+(sTSS+I1A-I1D)/2;j<(sTSS+I1A);j++){
      U5=U5+IN[i][j];
    }//for j
  }//for i
  U3=unsFLC;
  for(i=I1D+(sTSS+I1A-I1D)*5/6;i<L;i++){
    U3=U3+sP[i];
  }//for i
  for(i=I1D;i<I1D+(sTSS+I1A-I1D)/2;i++){
    for(j=I1D+(sTSS+I1A-I1D)*5/6;j<(sTSS+I1A);j++){
      U3=U3+IN[i][j];
    }//for j
  }//for i
  for(i=I1D+(sTSS+I1A-I1D)/2;i<I1D+(sTSS+I1A-I1D)*2/3;i++){
    for(j=i+(sTSS+I1A-I1D)/3;j<(sTSS+I1A);j++){
      U3=U3+IN[i][j];
    }//for j
  }//for i
  
  if (Us==0){//output conditional Loc|Us>0 distribution for comparison with model, but also keep info for when Us=0
    Loc=Loc-100;
  }//else

  output.open (fname.c_str(), ios::app);
  if(!output) {
    cout << "cannot open f file. \n";
  }//!foutput
  output << Vol << " " << sFLC  << " " << Us << " " << Um << " " << U5 << " " << U3 << " " << Loc << " " << unsFLC <<  " " << s1FLC << endl;
  output.close();

}//OutputToFile


void FaalFunctie(int* sP,int* sS1,int** IN,int sFLC,int unsFLC,int s1FLC,int STATE,double* p,double* p_group,int*** r,int* r_membernumber,int*sizecon,int** connect,int** reaction2group,int counter,double psum,int NR,int N_group,double t,int groupnumber,int currentreaction){
  
  int i,j;
  cout << endl << "counter= " << counter << " t=" << t << " groupnumber=" << groupnumber << " currentreaction=" << currentreaction << endl;
  cout << "sFLC= " << sFLC << endl <<  "unsFLC= " << unsFLC << endl << "s1FLC= " << s1FLC << endl << "STATE= " << STATE << endl;

  // cout << "*r " << endl;
  // for(i=0;i<NR;i++){
  //   if (r_membernumber[i]>0){
  //     cout << *r[i][0] << " ";
  //   }
  //   else{
  //     cout << "na ";
  //   }
  // }//for i
  // cout << endl;

  cout << "sP[] " << endl;
  for(i=0;i<L;i++){
    cout << sP[i] << " ";
  }//for i
  cout << endl;
  
  cout << "sS1[] " << endl;
  for(i=0;i<L;i++){
    cout << sS1[i] << " ";
  }//for i
  cout << endl;

  cout << "IN[] " << endl;

  // for(i=0;i<L;i++){
  //   for(j=0;j<L;i++){
  //     cout << IN[i][j] << " ";
  //   }//for j
  //   cout << endl;
  // }//for i
  // cout << endl;

   cout << "p[] "<< endl;
   //TEMP
   cout << p[2*L*L+2*L] << endl;
   //TEMP
  // for(i=0;i<NR;i++){
  //   cout << p[i] << " ";
  // }//for i
  // cout << endl;

  cout << "pgroup[] "<< endl;
  for(i=0;i<N_group;i++){
    cout << p_group[i] << " "; 
  }//for
  cout << endl << "psum = "  << psum << endl;
  for (i=0;i<sizecon[currentreaction];i++){
    cout << connect[currentreaction][i] << " ";
  }//for
  cout << endl << " r2g="<< reaction2group[currentreaction][0] << " "  << reaction2group[currentreaction][1] << " p[currentreaction]=" << p[currentreaction]<<  endl;
}//FaalFunctie
//temp


int main (int argc, char* argv[]) {
  //int main () {
  string date,sn,hname,fname,pname,sname,ename;
  ofstream output;
  int i,j,counter=0;
  int N_group,groupnumber,currentreaction;
  int sFLC=0,unsFLC=0,s1FLC=0,STATE=1;
  double t=0,tau,r1,pmin;
  double Vol,Nloci,beta,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9;

  //temp
  bool faal =0;
  //temp

  int NR = 2*L*L+2*L+7;//number of reactions
 
  date=string(argv[1]);
  sn=string(argv[2]); 
  Vol=atof(argv[3]); 

  //Model parameters
  Nloci=2.5;
  beta=31; //unit: pL^-1, exp determined
  k0 = 0.1;//atof(argv[4]);//elongation rate
  k1 = 0.0015;//atof(argv[5]);//splicing (intron processing) rate
  k2 = 0.067;//site s^-1 intronic RNA degradation 5'> 3': estimated
  k3 = 0.067;//site s^-1 intronic RNA degradation 3'> 5': estimated
  k4 = 0.02;//polyadenylation/P drop off, fixed parameter from literature Neugebauer PLoSB, Singer NatSMB
  k6 = 0.0005;//0.0003*Vol;//RNA release rate: manually fitted
  k7 = 0.000033;//FLC mRNA degradation rate t1/2=5.9h: exp determined
  k8 = 0.00005;//k_on: OFF>ON: manually fitted
  k9 = 0;//0.0003;//k_off: ON>OFF: manually fitted

  k5 = beta*Vol*k7/Nloci;//sense initiation rate
  
  cout << "FCA_52.cpp" << endl;
  
  pmin=SetPmin(k0,k1,k2,k3,k4,k5,k6,k7,k8,k9,N_group);// this also determines  N_group
  
  cout << "Time= " << Time << "s" << endl; 
  cout << "NR= " << NR << endl;
  cout << "pmin= " << pmin << endl;
  cout << "N_group= " << N_group << endl;
  cout << "L= " << L << " length of gene (sites), PolII footprint= 30bp" << endl;
  cout << "I1A (distance from sTSS, not site number!)= " << I1A << endl;
  cout << "I1D (site number)= " << I1D << endl;
  cout << "sTSS= " << sTSS << endl;
  cout << "asTSS= " << asTSS << endl;
  cout << "sTES= " << sTES << endl;
  cout << "cell Volume= " << Vol << endl;
  cout << "k0=elongation rate " << k0 <<  endl;//from Pol II fit, chromatinRNA fit
  cout << "k1=CoT splicing sense in1 rate " << k1 << endl;//estimated from smFISH
  cout << "k2=Intronic RNA degradation 5'> 3'(site/s) " << k2 << endl;
  cout << "k3=Intronic RNA degradation 3'> 5'(site/s) " << k3 << endl;
  cout << "k4=PolII drop off / pA s " << k4 << endl;//literature
  cout << "k5=sP s firing rate= " << k5 <<endl;//smFISH extracted
  cout << "k6=(s) RNA export rate " << k6 << endl;//estimated from fit smFISH
  cout << "k7=FLC RNA degradation rate " << k7 << endl;//smFISH measured
  cout << "k8=OFF>ON STATE rate " << k8 << endl;
  cout << "k9=ON>OFF STATE rate " << k9 << endl;
  
  date.append("_"); 
  hname = date;
  fname = date;
  pname = date;
  hname.append("h");
  fname.append("f");
  pname.append("p");
  hname.append(sn);
  fname.append(sn);
  pname.append(sn);
  hname.append(".txt");
  fname.append(".txt");
  pname.append(".txt");

  //step(0) initialize reactants and other arrays
  int * sP=0;
  sP = new int [L];//number of sense PolII at sites
  int * sS1=0;
  sS1 = new int [L];//number of sense PolII at sites
  int ** IN=0;
  IN = new int* [L];//number of intron lariats with end at sites i,j
  for (i=0;i<L;i++){
    IN[i]= new int [L];
    sP[i]=0;
    sS1[i]=0;
  }//for i
  for (i=0;i<L;i++){
    for (j=0;j<L;j++){
      IN[i][j]=0;
    }//for j
  }//for i

  //group2reaction has indices of propensities as entries ordered by the group they belong t
  int ** group2reaction = 0;
  group2reaction = new int * [N_group];
  for (i=0;i<N_group;i++){
    group2reaction [i]= new int [NR];
    for (j=0;j<NR;j++){
      group2reaction[i][j]=-1;
    }//for j
  }//for i

  //entries are number of reactions in particular group
  int * g_membernumber=0;
  g_membernumber = new int [N_group];
  for(i=0;i<N_group;i++){
    g_membernumber[i]=0;//initially
  }//for i

  int* r_membernumber=0;
  r_membernumber = new int [NR];
  SetR_membernumber (r_membernumber,L);

  int*** r = 0;
  r = new int** [NR];
  for(i=0;i<NR;i++){
    if (r_membernumber[i]>0){
      r[i]= new int* [r_membernumber[i]];
    }//if
    else{
      r[i] = 0;
    }//else
  }//for i
 
  SetR(r,L,sP,sS1,IN,&sFLC,&unsFLC,&s1FLC,&STATE,sTSS,sTES);
  char* reactiontype = 0;
  reactiontype = new char [NR];
  SetRT(reactiontype,L,I1A,I1D,sTSS,sTES);

  double** k = 0;
  k = new double* [NR];
  for(i=0;i<NR;i++){
    k[i]= new double [1];
    k[i][0]=0;
  }//for i
 
  SetK(k,L,k0,k1,k2,k3,k4,k5,k6,k7,k8,k9);

  double* p = 0;
  p = new double [NR];
 
  int** reaction2group = 0;
  reaction2group = new int * [NR];
  for (i=0;i<NR;i++){
    reaction2group[i]= new int [2];//[0]=group, [1]=position
    reaction2group[i][0]=-1;//initially
    reaction2group[i][1]=-1;
  }//for i
 
  SetP(p,k,r,reactiontype,g_membernumber,group2reaction,NR, pmin,reaction2group);

  int * sizecon = 0;
  sizecon = new int [NR];
  SetSizecon(sizecon,L,sTSS,sTES);

  int ** connect =0;
  connect = new int * [NR];
  for(i=0;i<NR;i++){
    connect[i] = new int [sizecon[i]];
  }// i
  SetCon(connect,L,sTSS,sTES,I1D,I1A);

  double * p_group = 0;
  p_group = new double[N_group];
  SetPgroup (p,p_group, N_group,g_membernumber,group2reaction);
  
  double psum;
  psum=SetPsum(p_group,N_group); 

  do {
    
    //  cout << "loop0" << endl;
    r1=mtrand1.randDblExc();//step 1

    tau=(1/psum)*log(1/r1);//step 2 
 
    do{
      groupnumber = GetGroupNumber(p_group,psum,N_group, faal);//step 3a 

      //quality control
      if((groupnumber < 0) || (groupnumber >= N_group) || g_membernumber[groupnumber]==0 || p_group[groupnumber]==0  ){
	//faal =1;
	cout << "groupnumber assignment gaat fout groupnumber=" << groupnumber <<  endl;
      }//if
      //quality control
    }while((g_membernumber[groupnumber]==0) || (p_group[groupnumber]==0));

    //   cout << "loop2" << endl;

    currentreaction = RejectionSample(group2reaction,p,g_membernumber[groupnumber],groupnumber,pmin,faal,t);//step 3b  
    //cout << currentreaction << endl;

    //quality control
    if(currentreaction < 0 || currentreaction > NR){
      faal =1;
      cout << "current reaction gaat fout reactionindex=" << currentreaction <<  endl;
    }//if
    //quality control

    PerformReaction(currentreaction,r,connect,sizecon,sP,sS1,IN,sFLC,unsFLC,s1FLC,STATE,t,faal);//step 4 

    UpdateP1(currentreaction,p,k,r,reactiontype,connect,sizecon,group2reaction,g_membernumber,reaction2group,p_group,N_group,pmin,psum);//step 5 and 6a
    //    cout << "loop3" << endl;

    //quality control
    if ( (psum<=0) ){
      faal = 1;
      cout << "psum gaat fout pmin*pow(2.0,N_group)= " << pmin*pow(2.0,N_group);
    }//quality control

    //temp
    if (faal == 1){
      FaalFunctie(sP,sS1,IN,sFLC,unsFLC,s1FLC,STATE,p,p_group,r,r_membernumber,sizecon,connect,reaction2group,counter,psum,NR,N_group,t,groupnumber,currentreaction);

      //faal=0;
      break;
    }//if
    //temp

    t=t+tau;
   
    counter++;

  }while(t<Time); 
  
  OutputToFile(sP,sS1,IN,sFLC,unsFLC,s1FLC,Vol,fname);
  
  return 0;

}//main


