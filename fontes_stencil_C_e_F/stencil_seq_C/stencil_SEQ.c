#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>


//cc stencil.c -o stencil_pgi_bro.x -lm
//                             $0          n   energy niters
//time aprun -b -n 1 -d 1 -N 1 stencil_pgi 100 1      200 
//Xmitswan.sh 1 1 01 BROADWELL MPI


// row-major order
#define ind(i,j) (j)*(n+2)+i



int main(int argc, char **argv) 
{
  int n = atoi(argv[1]); // nxn grid :: 100x100=10000
  int energy = atoi(argv[2]); // energy to be injected per iteration :: 1
  int niters = atoi(argv[3]); // number of iterations :: 200
  double *aold = (double*)calloc(1,(n+2)*(n+2)*sizeof(double)); // 1-wide halo zones!
  double *anew = (double*)calloc(1,(n+2)*(n+2)*sizeof(double)); // 1-wide halo-zones!
  double t=0.0, tinit=0.0, tstencil=0.0, tenergy=0.0, tresto=0.0;
  int i, iter, j;

  MPI_Init(NULL, NULL);
  t=-MPI_Wtime();

  tinit-=MPI_Wtime();
  for(i=0; i<(n+2)*(n+2); i++)
  {
    aold[i]=0.0;
    anew[i]=0.0;
  }

  #define nsources 3
  int sources[nsources][2] = {{n/2,n/2}, {n/3,n/3}, {n*4/5,n*4/5}};
  tinit+=MPI_Wtime();


  for(iter=0; iter<niters; iter+=2) 
  {
    
    tstencil-=MPI_Wtime();
    for(j=1; j<n+1; ++j) 
      for(i=1; i<n+1; ++i) 
        anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
    tstencil+=MPI_Wtime();

    tenergy-=MPI_Wtime();
    for(i=0; i<nsources; ++i) 
      anew[ind(sources[i][0],sources[i][1])] += energy; // heat source
    tenergy+=MPI_Wtime();

    tstencil-=MPI_Wtime();
    for(j=1; j<n+1; ++j) 
      for(i=1; i<n+1; ++i) 
        aold[ind(i,j)] = anew[ind(i,j)]/2.0 + (anew[ind(i-1,j)] + anew[ind(i+1,j)] + anew[ind(i,j-1)] + anew[ind(i,j+1)])/4.0/2.0;

    tstencil+=MPI_Wtime();
    

    tenergy-=MPI_Wtime();
    for(i=0; i<nsources; ++i) 
      aold[ind(sources[i][0],sources[i][1])] += energy; // heat source
    tenergy+=MPI_Wtime();
  }

  tenergy-=MPI_Wtime();
  double heat=0.0; // total heat in system
  heat=0.0;
  for(j=1; j<n+1; ++j) 
    for(i=1; i<n+1; ++i) 
      heat += aold[ind(i,j)];
  tenergy+=MPI_Wtime();
  t+=MPI_Wtime();
  tresto=t-(tstencil+tinit+tenergy);

  printf("SEQ_C       execucao com %2d procs  heat= %f  Tempos: Total= %f  Stencil= %f  Troca= %f  Energy= %f  Init= %f  Reduce= %f  Resto= %f\n",
     1, heat, t, tstencil, 0.0, tenergy, tinit, 0.0, tresto);
  MPI_Finalize();
  
}
