#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#define m 1000
#define n 1000
#define T 2
#define dt 0.01
#define dx 0.1
#define D 0.1
//=========================
void DisplayMatrix(float *A, int row, int col)
{
int i,j;
for(i=0;i<row;i++){
for(j=0;j<col;j++) printf(" %f",*(A+i*col+j));
printf("\n");
}
}
//=========================
void Write2File(float *C)
{
FILE *result=fopen("result2.txt", "a");
int i,j;
for(i=0;i<m;i++)
{
for(j=0;j<n;j++)
{
fprintf(result, "%lf\t", *(C+i*n+j));
}
fprintf(result, "\n");
}
fclose(result);
}
//=========================
void KhoiTao(float *C)
{
int i,j;
for ( i = 0 ; i < m ; i++ )
for ( j = 0 ; j < n ; j++ ){
if (i>=(m/2-5)&&i<(m/2+5)&&j>=(n/2-5)&&j<(n/2+5))
*(C+i*n+j) = 80.0;
else
*(C+i*n+j) = 25.0;
}
}
//=========================
void FD(float *C, float *dC,int Mc) {
int i, j;
float c,u,d,l,r;
for ( i = 0 ; i < Mc ; i++ )
for ( j = 0 ; j < n ; j++ )
{
c = *(C+i*n+j);
u = (i==0) ? *(C+i*n+j) : *(C+(i-1)*n+j);
d = (i==m-1) ? *(C+i*n+j) : *(C+(i+1)*n+j);
l = (j==0) ? *(C+i*n+j) : *(C+i*n+j-1);
r = (j==n-1) ? *(C+i*n+j) : *(C+i*n+j+1);
*(dC+i*n+j) = (1/(dx*dx))*(u+d+l+r-4*c);
}
}
//=========================
int main(int argc,char **argv)
{
//khaibaobien
int i,j, count;
float t;
time_t t1,t2;
int NP,rank,Mc;
MPI_Init(&argc,&argv );
MPI_Comm_size(MPI_COMM_WORLD,&NP);
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Status thongbao;
float *Cs,*dCs,Cu,Cd,Cl,Cr,*dC,*C;
t =0;
Mc = m/NP;
dCs =(float *)malloc ((Mc*n)*sizeof(float));
C = (float *) malloc ((m*n)*sizeof(float));
dC = (float *) malloc ((m*n)*sizeof(float));
Cs = (float *) malloc ((m*n)*sizeof(float));

//phantandulieu
if(rank ==0)
  for(i=0;i<m;i++){
  	for ( j = 0 ; j < n ; j++ )
  		*(C+i)=25.0;
  		 MPI_Scatter (C, Mc*n, MPI_FLOAT,Cs, Mc*n, MPI_FLOAT, 0,MPI_COMM_WORLD);
  }
//khoitao
KhoiTao(C);
Write2File(C);
//printf("Gia tri khoi tao:\n");
//DisplayMatrix(C, m, n);
count=0;
t1=time(NULL);
//tinhtoan
while (t<=T)
{
if(rank==(NP-1)){
    Cu = 80.0;
    MPI_Send(Cs+Mc-1,1,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD);
  }
  else if(rank!= (NP-1)){
    MPI_Recv(&Cu,1,MPI_FLOAT,rank+1,rank, MPI_COMM_WORLD,&thongbao);
    MPI_Send(Cs,1,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD);
  }
  else{
    MPI_Recv(&Cu,1,MPI_FLOAT,rank+1,rank, MPI_COMM_WORLD,&thongbao);
    
  }
  if(rank==(NP-1)){
    Cd = 25.0;
    MPI_Send(Cs+Mc-1,1,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD);
  }
  else if(rank !=NP-1){
    MPI_Recv(&Cd,1,MPI_FLOAT,rank+1,rank, MPI_COMM_WORLD,&thongbao);
    MPI_Send(Cs,1,MPI_FLOAT,rank-1,rank-1,MPI_COMM_WORLD);
  }
  else{
    MPI_Recv(&Cd,1,MPI_FLOAT,rank+1,rank, MPI_COMM_WORLD,&thongbao);
  }
FD(C, dC,Mc);
for ( i = 0 ; i < m ; i++ )
for ( j = 0 ; j < n ; j++ )
*(C+i*n+j) = *(C+i*n+j) + dt*(*(dC+i*n+j));
t=t+dt;
count=count+1;
if (count%5==0) Write2File(C);
}
t2=time(NULL);
//FD(C, dC);
//printf("Gia tri cuoi cung:\n");
//DisplayMatrix(C, m, n);
//
printf ("\tThe Calculation time:%ld\n",(long)(t2 - t1));
MPI_Gather ( Cs, Mc*n, MPI_FLOAT,C, Mc*n, MPI_FLOAT, NP-1,MPI_COMM_WORLD);
  // in ra KQ
if(rank==NP-1)
  for(i=0;i<m;i++) printf("%f\n",*(C+i) );
  MPI_Finalize();
return 0;
}