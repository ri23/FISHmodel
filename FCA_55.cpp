//Robert Ietswaart
//20170223
//FCA_55.cpp

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
#include "FCA_55_table.h"

using namespace std;

const double Time = 10*24*60*60;//Time in s
const int cap =1000;//artificial cap that FLC levels shouldn't exceed

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

double SetPmin (double k0,double k1,double k2, double k3,int& Ng){
 
  int arraysize=4;//arraysize is number of rates involved
  double kmin, k_sort[arraysize];
  k_sort[0]=k0;
  k_sort[1]=k1;
  k_sort[2]=k2;
  k_sort[3]=k3;//only when k3>0
  sort(k_sort,k_sort+arraysize);
  kmin= k_sort[0];

  k_sort[0]=k0;
  k_sort[1]=k1*cap;
  k_sort[2]=k2;
  k_sort[3]=k3;//only when k3>0

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
 
double ProductionP(double k, int** r){
  double pnew;
  if(*r[0]==1){
    pnew=k;
  }//if new site not occupied by P2
  else{
    pnew=0;
  }//if new site occupied, can't fire
  //r[0]= STATE
  return pnew;
}//ProductionP

//site is reaction index
double CalculateP (double* k, int*** r, char RT, int site){
  double pnew;
  if (RT==0){
    pnew=ProductionP(k[site],r[site]);
  }//zeroth order reaction (A production)
  else if (RT==1){
    pnew = k[site]*(*r[site][0]);//first order reaction   
  }//first order reaction (B degradation) / F switch to OFF
  else if (RT==2){    
    pnew = k[site]*(1-*r[site][0]);
  }// if D switch to ON STATE
  else if (RT==3){    
    pnew = k[site];
  }// if C void>void
  else{
    cout << "faal, dit reactietype is niet toegestaan " << RT << " " << k[site] << " site="<< site << endl;
    pnew=0;//temp faal oplossing, deze shit mag eenvoudig niet voorkomen
  }// else
  // cout << pnew << endl;
  return pnew;
}//double CalculateP


void SetP (double * p, double* k, int*** r, char* RT,  int * g_membernumber,int ** group2reaction, int NR ,double pmin, int** reaction2group){
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
  // if (pg[group]<((pmin*8)/10)){//Shabby solution for problem of floating point values
  //   pg[group]=0;
  // }//if
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
void PerformReaction(int m, int*** r, int** connect, int* sizecon,int& sFLC, int& STATE,bool& fout){
  int i;//positions on grid
  if (m==0){
    sFLC++;
    if (sFLC>cap){
      cout << endl << "cA " << "sFLC= " << sFLC << " m=" << m <<   endl; 
      fout=1;
    }//if
  }//if m==0
  else if (m==1){
    sFLC--;
    if (sFLC<0){
      cout << endl << "cB " << "sFLC= " << sFLC << " m=" << m <<   endl; 
      fout=1;
    }//if
  }//else if m==1
  else if (m==2){
    //token reaction to keep program running
  }//else if m==2
  //temp
  else if (m == 3){
    STATE=1;
  }//D 2 1: OFF > ON
  else if (m == 4){
    STATE=0;
  }//E 3 1: ON > OFF
  else{
    cout << "faal Perf reacton m=" << m << endl;
  }//else
  //temp
}//PerformReaction

//you only get in here if propensity actually changes
void UpdateP1 (int m, double* p, double* k, int*** r, char* rt, int** connect , int* sizecon,int** group2reaction, int* g_membernumber, int** reaction2group, double* pg, int Ng, double pmin, double& ps){
  int i, site;
  double pnew;
  //cout << "m= " << m << endl;
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


void FaalFunctie(int sFLC,double* p,double* p_group,int*** r,int* r_membernumber,int*sizecon,int** connect,int** reaction2group,int counter,double psum,int NR,int N_group,double t,int groupnumber,int currentreaction){
  
  int i;
  cout << endl << "counter= " << counter << " t=" << t << " groupnumber=" << groupnumber << " currentreaction=" << currentreaction << endl;
  cout <<  "sFLC= " << sFLC << endl;

  cout << "*r " << endl;
  for(i=0;i<NR;i++){
    if (r_membernumber[i]>0){
      cout << *r[i][0] << " ";
    }
    else{
      cout << "na ";
    }
  }//for i
  cout << endl;

  cout << "p[] "<< endl;
  for(i=0;i<NR;i++){
    cout << p[i] << " ";
  }//for i
  cout << endl;

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
  string date,sn,hname;
  ofstream output;
  int i,j,counter=0;
  int N_group,groupnumber,currentreaction;
  int sFLC=0,STATE=1;
  double t=0,tau,r1,pmin;
  double beta,Vol,Vol_av,k0,k1,k2,k3,bs,N_loci;

  int NR = 5;//number of reactions

  //temp
  bool faal =0;
  //temp

  date=string(argv[1]);
  sn=string(argv[2]);
  Vol=atof(argv[3]); 
  N_loci=atof(argv[4]);
  beta=31;
  k1 = 0.000033;//mRNA degradation rate
  k3 = 0.1;//k_off: ON>OFF: manually set to ensure koff>>kon
  //SCENARIO1: burst size scales with volume, not burst frequency (k_on)
  bs=7;//burst size (for average volume): set to vary manually
  Vol_av=1.8;//average volume (pL)
  k0 = bs*Vol/Vol_av*k3;//production rate in ON state: p_on~Vol
  k2 = beta*Vol_av*k1*k3/(k3*bs*N_loci-beta*Vol*k1);//k_on: OFF>ON: independent of vol
  bs = k0/k3;//burst size for this particular cell/Vol
  //SCENARIO2: burst frequency (k_on) scales with volume, not burst size
  //bs= 6;//burst size: set to vary manually
  //k0= bs*k3;//production rate in ON state: p_on
  //k2= beta*Vol*k1*k3/(k3*bs*N_loci-beta*Vol*k1);//k_on dependent of vol and N_loci
  //SCENARIO3: Poisson production
  //don't forget to adjust in SetPmin to array size 3 and exclude k3(=0) to prevent segmentation fault
  //bs = 1;//does not have any meaning here
  //k0 = beta*Vol*k1/N_loci;//p_on: only relevant parameter here
  //k2 = 1;//k_on: OFF>ON: does not have any meaning here as systems never reaches OFF state
  //k3 = 0;//k_off: ON>OFF: set to zero

  pmin=SetPmin(k0,k1,k2,k3,N_group);// this also determines  N_group

  cout << "FCA_55.cpp" << endl;
  cout << "Time= " << Time << "s" << endl; 
  cout << "NR= " << NR << endl;
  cout << "pmin= " << pmin << endl;
  cout << "N_group= " << N_group << endl;
  cout << "cell Volume= " << Vol << endl;
  cout << "N_loci= " << N_loci << endl;
  cout << "beta= " << beta << endl;
  cout << "burst size= " << bs << endl;
  cout << "k0=FLC mRNA production rate per locus= " << k0 <<  endl;
  cout << "k1=FLC mRNA degradation rate " << k1 << endl;
  cout << "k2=OFF>ON STATE rate=burst freq " << k2 << endl;
  cout << "k3=ON>OFF STATE rate " << k3 << endl;
  cout << "koff/kon ratio= " << k3/k2 << endl;

  date.append("_"); 
  hname = date;
  hname.append("f");
  hname.append(sn);
  date.append(".txt");
  hname.append(".txt");

  //step(0) initialize reactants and other arrays
 
  //group2reaction has indices of propensities as entries ordered by the group they belong to
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
  SetR_membernumber (r_membernumber);

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
 
  SetR(r,&sFLC,&STATE);
  char* reactiontype = 0;
  reactiontype = new char [NR];
  SetRT(reactiontype);

  double* k = 0;
  k = new double [NR];
  for(i=0;i<NR;i++){
    k[i]=0;
  }//for i
 
  SetK(k,k0,k1,k2,k3);

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
  SetSizecon(sizecon);

  int ** connect =0;
  connect = new int * [NR];
  for(i=0;i<NR;i++){
    connect[i] = new int [sizecon[i]];
  }// i
  SetCon(connect);

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
    }while((g_membernumber[groupnumber]==0) && (p_group[groupnumber]==0));

    //   cout << "loop2" << endl;

    currentreaction = RejectionSample(group2reaction,p,g_membernumber[groupnumber],groupnumber,pmin,faal,t);//step 3b  
    //cout << currentreaction << endl;

    //quality control
    if(currentreaction < 0 || currentreaction > NR){
      faal =1;
      cout << "current reaction gaat fout reactionindex=" << currentreaction <<  endl;
    }//if
    //quality control

    PerformReaction(currentreaction,r,connect,sizecon,sFLC,STATE,faal);//step 4 

    UpdateP1(currentreaction,p,k,r,reactiontype,connect,sizecon,group2reaction,g_membernumber,reaction2group,p_group,N_group,pmin,psum);//step 5 and 6a
    //    cout << "loop3" << endl;

    //quality control
    if ( (psum<=0) ){
      faal = 1;
      cout << "psum gaat fout pmin*pow(2.0,N_group)= " << pmin*pow(2.0,N_group);
    }//quality control

    //temp
    if (faal == 1){
      FaalFunctie(sFLC,p,p_group,r,r_membernumber,sizecon,connect,reaction2group,counter,psum,NR,N_group,t,groupnumber,currentreaction);
      //faal=0;
      break;
    }//if
    //temp

    t=t+tau;
   
    counter++;

  }while(t<Time); 
  

  output.open (hname.c_str(), ios::app);
  if(!output) {
    cout << "cannot open file. \n";
  }//if
  output << Vol << " " << sFLC << endl;
  output.close();
 
  return 0;

}//main


