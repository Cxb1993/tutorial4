#include "sor.h"
#include <math.h>
#include "parallel.h"


void sor(
         double omg,
         double dx,
         double dy,
         int    il,
         int    ir,
         int 	jb,
         int	jt,
         int rank_l,
         int rank_r,
         int rank_b,
         int rank_t,
         double *bufSend,
         double *bufRecv,
         MPI_Status *status,
         int chunk,
         double **P,
         double **RS,
         int myrank,
         int imax,
         int jmax,
         double *res,
         int omg_i,
         int omg_j,
         int iproc,
         int jproc
         ) {
    
    int i,j;
    double rloc;
    double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
    
    double glob_res;
    
    
    
    
    /* SOR iteration */
    
    
    
    for(j = jb+1; j <= jt; j++) {
        for(i = il+1; i<=ir; i++) {
            P[i][j] = (1.0-omg)*P[i][j]
            + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
        }
    }
    
    /*Boundary values*/
    
    if(omg_i==1){
        for (j = jb; j <= jt+1; j++){
            /*U values around left boundary*/
            P[il][j] = P[il+1][j];
        }
    }
    if(omg_i==iproc){
        for (j = jb; j <= jt+1; j++){
            /*U values around left boundary*/
            P[ir+1][j] = P[ir][j];
        }
    }
    
    if(omg_j==1){
        for (i = il; i <= ir+1; i++){
            /*U values around left boundary*/
            P[i][jb] = P[i][jb+1];
        }
    }
    if(omg_j==jproc){
        for (i = il; i <= ir+1; i++){
            /*U values around left boundary*/
            P[i][jt+1] = P[i][jt];
        }
    }
    
    /*Passing the pressure values*/
    pressure_comm(P,il,ir,jb,jt ,rank_l,rank_r,rank_b,rank_t,bufSend,bufRecv,status,chunk);
        
    /* compute the residual */
    rloc = 0;
    for(j = jb+1; j <= jt; j++) {
        for(i = il+1; i <= ir; i++) {
            rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
            ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
        }
    }
    MPI_Reduce(&rloc, &glob_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if ( myrank == 0 ) {
        
        glob_res = glob_res/(jmax*imax) ;
        glob_res = sqrt(glob_res) ;
        /*printf("my rank = %i, global = %f\n", myrank, glob_res);*/
    }
    MPI_Bcast (&glob_res,1,MPI_DOUBLE,0,MPI_COMM_WORLD) ;
    *res=glob_res;
    
}
