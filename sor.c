#include "sor.h"
#include <math.h>
#include "parallel.h"


void sor(
         double omg,
         double dx,
         double dy,
         int    ir,
         int    il,
         int 	jt,
         int	jb,
         double **P,
         double **RS,
         int myrank,
         int rank_l,
         int rank_r,
         int rank_b,
         int rank_t,
         double *bufSend,
         double *bufRecv,
         MPI_Status *status,
         int chunk,
         int imax,
         int jmax,
         
         ) {
    
    int i,j;
    double rloc;
    double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));
    double glob_res  ;
    
    
    /*Boundary values*/
    
    for (j = jb; j <= jt; j++){
        
        if ( (il-1) == 0 ){
            /*U values around left boundary*/
            P[il-1][j] = P[il][j];
        }
        else if( ir == imax){
            /*U values around right boundary*/
            P[ir+1][j] = P[ir][j];
        }
    }
    
    for (i = il; i <= ir; i++){
        
        if ( (jb-1) == 0 ){
            /*P values around bottom boundary*/
            P[i][jb-1] = P[i][jb];
        }
        else if( jt == jmax){
            /*P values around top boundary*/
            P[i][jt+1] = P[i][jt];
        }
    }
    
    
    
    
    /* SOR iteration */
    
    
    
    for(j = jb; j <= jt; j++) {
        for(i = il; i<=ir; i++) {
            P[i][j] = (1.0-omg)*P[i][j]
            + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
        }
    }
    
    
    /*Passing the pressure values*/
    pressure_comm(P,il,ir,jb,jt ,rank_l,rank_r,rank_b,rank_t,bufSend,bufRecv,status,chunk );
    
    
    
    
    
    /* compute the residual */
    rloc = 0;
    for(j = jb; j <= jt; j++) {
        for(i = il; i <= ir; i++) {
            rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
            ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
        }
    }
    
    MPI_Reduce(&rloc, &glob_res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    
    if ( myrank == 0 ) {
        
        glob_res = glob_res/(jmax*imax) ;
        glob_res = sqrt(glob_res) ;
        
        MPI_Bcast (&glob_res,1,MPI_DOUBLE,0,MPI_COMM_WORLD) ;
        
        
    }
    
    
}
