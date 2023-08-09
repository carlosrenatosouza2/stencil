/*
 * Copyright (c) 2012 Torsten Hoefler. All rights reserved.
 *
 * Author(s): Torsten Hoefler <htor@illinois.edu>
 *
 */

#include "stencil_par.h"

int main(int argc, char **argv) {
  double t=0.0, tinit=0.0, tstencil=0.0, tenergy=0.0;
  double ttroca=0.0, tred=0.0, tresto=0.0;
  MPI_Init(&argc, &argv); 
  int r,p;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &r);
  MPI_Comm_size(comm, &p);
  int n, energy, niters, i, j, iter;



  if (r==0) {
      // argument checking
      if(argc < 4) {
          if(!r) printf("usage: stencil_mpi <n> <energy> <niters>\n");
          MPI_Finalize();
          exit(1);
      }

      n = atoi(argv[1]); // nxn grid
      energy = atoi(argv[2]); // energy to be injected per iteration
      niters = atoi(argv[3]); // number of iterations
      
      // distribute arguments
      int args[3] = {n, energy, niters};
      MPI_Bcast(args, 3, MPI_INT, 0, comm);
  }
  else {
      int args[3];
      MPI_Bcast(args, 3, MPI_INT, 0, comm);
      n=args[0]; energy=args[1]; niters=args[2];
  }

  MPI_Comm shmcomm;
  MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);

  int sr, sp; // rank and size in shmem comm
  MPI_Comm_size(shmcomm, &sp);
  MPI_Comm_rank(shmcomm, &sr);

  // this code works only on comm world!
  if(sp != p) MPI_Abort(comm, 1);

  int pdims[2]={0,0};
  MPI_Dims_create(sp, 2, pdims);
  int px = pdims[0];
  int py = pdims[1];
  //if(!r) printf("processor grid: %i x %i\n", px, py);
  t-=MPI_Wtime(); // take time
  tinit-=MPI_Wtime(); // take time

  // determine my coordinates (x,y) -- r=x*a+y in the 2d processor array
  int rx = r % px;
  int ry = r / px;
  // determine my four neighbors
  int north = (ry-1)*px+rx; if(ry-1 < 0)   north = MPI_PROC_NULL;
  int south = (ry+1)*px+rx; if(ry+1 >= py) south = MPI_PROC_NULL;
  int west= ry*px+rx-1;     if(rx-1 < 0)   west = MPI_PROC_NULL;
  int east = ry*px+rx+1;    if(rx+1 >= px) east = MPI_PROC_NULL;
  // decompose the domain
  int bx = n/px; // block size in x
  int by = n/py; // block size in y
  int offx = rx*bx; // offset in x
  int offy = ry*by; // offset in y

  int size = (bx+2)*(by+2); // process-local grid (including halos (thus +2))
  double *mem;
  MPI_Win win;
  MPI_Win_allocate_shared(2*size*sizeof(double), 1, MPI_INFO_NULL, shmcomm, &mem, &win);


  double *tmp;
  double *anew=mem; // each rank's offset
  double *aold=mem+size; // second half is aold!

  double *northptr, *southptr, *eastptr, *westptr;
  double *northptr2, *southptr2, *eastptr2, *westptr2;
  MPI_Aint sz;
  int dsp_unit;
  if(north != MPI_PROC_NULL) 
  {
    MPI_Win_shared_query(win, north, &sz, &dsp_unit, &northptr);
    northptr2 = northptr+size;
  }

  if(south != MPI_PROC_NULL) 
  {
    MPI_Win_shared_query(win, south, &sz, &dsp_unit, &southptr);
    southptr2 = southptr+size;
  }

  if(east != MPI_PROC_NULL) 
  {
    MPI_Win_shared_query(win, east, &sz, &dsp_unit, &eastptr);
    eastptr2 = eastptr+size;
  }

  if(west != MPI_PROC_NULL) 
  {
    MPI_Win_shared_query(win, west, &sz, &dsp_unit, &westptr);
    westptr2 = westptr+size;
  }  
  // initialize three heat sources
  #define nsources 3
  int sources[nsources][2] = {{n/2,n/2}, {n/3,n/3}, {n*4/5,n*4/5}};
  int locnsources=0; // number of sources in my area
  int locsources[nsources][2]; // sources local to my rank
  for (i=0; i<nsources; ++i) // determine which sources are in my patch
  { 
    int locx = sources[i][0] - offx;
    int locy = sources[i][1] - offy;
    if(locx >= 0 && locx < bx && locy >= 0 && locy < by) 
    {
      locsources[locnsources][0] = locx+1; // offset by halo zone
      locsources[locnsources][1] = locy+1; // offset by halo zone
      locnsources++;
    }
  }
  
  for(i=0; i<(bx+2)*(by+2); i++)
  {
    aold[i]=0.0;
    anew[i]=0.0;
  }

  MPI_Win_lock_all(0, win);
  tinit+=MPI_Wtime(); // take time

  for(iter=0; iter<niters; iter+=2) 
  {

    // anew <- stencil(aold)

    // exchange data with neighbors
    ttroca-=MPI_Wtime(); // take time
    MPI_Win_sync(win);
    MPI_Barrier(shmcomm);
    if(north != MPI_PROC_NULL)
      for(i=0; i<bx; ++i) aold[ind(i+1,0)] = northptr2[ind(i+1,by)]; // pack loop - last valid region
    if(south != MPI_PROC_NULL)
      for(i=0; i<bx; ++i) aold[ind(i+1,by+1)] = southptr2[ind(i+1,1)]; // pack loop
    if(east != MPI_PROC_NULL)
      for(i=0; i<by; ++i) aold[ind(bx+1,i+1)] = eastptr2[ind(1,i+1)]; // pack loop
    if(west != MPI_PROC_NULL)
      for(i=0; i<by; ++i) aold[ind(0,i+1)] = westptr2[ind(bx,i+1)]; // pack loop
    ttroca+=MPI_Wtime();

    // update grid points
    tstencil-=MPI_Wtime();
    for(j=1; j<by+1; ++j)
      for(i=1; i<bx+1; ++i)
        anew[ind(i,j)] = aold[ind(i,j)]/2.0 + (aold[ind(i-1,j)] + aold[ind(i+1,j)] + aold[ind(i,j-1)] + aold[ind(i,j+1)])/4.0/2.0;
    tstencil+=MPI_Wtime();
    
    // refresh heat sources
    tenergy-=MPI_Wtime();
    for(i=0; i<locnsources; ++i)
      anew[ind(locsources[i][0],locsources[i][1])] += energy; // heat source
    tenergy+=MPI_Wtime();

    // aold <- stencil(anew)

    ttroca-=MPI_Wtime(); // take time
    // exchange data with neighbors
    MPI_Win_sync(win);
    MPI_Barrier(shmcomm);
    if(north != MPI_PROC_NULL)
      for(i=0; i<bx; ++i) anew[ind(i+1,0)] = northptr[ind(i+1,by)]; // pack loop - last valid region
    if(south != MPI_PROC_NULL)
      for(i=0; i<bx; ++i) anew[ind(i+1,by+1)] = southptr[ind(i+1,1)]; // pack loop
    if(east != MPI_PROC_NULL)
      for(i=0; i<by; ++i) anew[ind(bx+1,i+1)] = eastptr[ind(1,i+1)]; // pack loop
    if(west != MPI_PROC_NULL)
      for(i=0; i<by; ++i) anew[ind(0,i+1)] = westptr[ind(bx,i+1)]; // pack loop
    ttroca+=MPI_Wtime();

    // update grid points
    tstencil-=MPI_Wtime();
    for(j=1; j<by+1; ++j)
      for(i=1; i<bx+1; ++i)
        aold[ind(i,j)] = anew[ind(i,j)]/2.0 + (anew[ind(i-1,j)] + anew[ind(i+1,j)] + anew[ind(i,j-1)] + anew[ind(i,j+1)])/4.0/2.0;
    tstencil+=MPI_Wtime();
    
    // refresh heat sources
    tenergy-=MPI_Wtime();
    for(i=0; i<locnsources; ++i)
      aold[ind(locsources[i][0],locsources[i][1])] += energy; // heat source
    tenergy+=MPI_Wtime();
  }

  tenergy-=MPI_Wtime();
  MPI_Win_unlock_all(win);
  double heat = 0.0;
  for(j=1; j<by+1; ++j)
    for(i=1; i<bx+1; ++i)
      heat += aold[ind(i,j)];
  tenergy+=MPI_Wtime();


  t+=MPI_Wtime();

  // get final heat in the system
  tred-=MPI_Wtime();
  double rheat;
  MPI_Allreduce(&heat, &rheat, 1, MPI_DOUBLE, MPI_SUM, comm);
  tred+=MPI_Wtime();
  tresto=t-(tstencil+ttroca+tenergy+tinit);

  if(!r)
    printf("SHM_C       execucao com %2d procs  heat= %f  Tempos: Total= %f  Stencil= %f  Troca= %f  Energy= %f  Init= %f  Reduce= %f  Resto= %f\n", 
	   p, rheat, t, tstencil, ttroca, tenergy, tinit, tred, tresto);

  MPI_Win_free(&win);
  MPI_Comm_free(&shmcomm);
  MPI_Finalize();
}
