#include "boundary_val.h"

/* ----------------------------------------------------------------------- */
/*                             Set Boundary Conditions                     */
/* ----------------------------------------------------------------------- */

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
                    int il,
                    int ir,
                    int jb,
                    int jt,
                    int imax,
                    int jmax,
                    double **U,
                    double **V
                    ) {
    
    int i,j;
    
    /*Set boundary values along the columns*/
    for (j = jb-1; j <= jt+1; j++){
        
        if ( (il) == 0 ){
            /*U velocities on left boundary*/
            U[il-1][j] = 0;
            /*V velocities on left boundary (interpolated value)*/
            V[il-1][j]=-1*V[il][j];
        }
        if( ir + 1 == imax){
            /*U velocities on right boundary*/
            U[ir][j] = 0;
            /*V velocities on right boundary (interpolated value)*/
            V[ir+1][j]=-1*V[ir][j];
        }
    }
    
    /*Set boundary values along the rows*/
    for (i = il-1; i <= ir+1; i++){

        if ( (jb) == 0 ){
            /*V velocities on bottom boundary*/
            V[i][jb-1] = 0;
            /*U velocities on bottom boundary (interpolated value)*/
            U[i][jb-1]=-1*U[i][jb];
            
        }
        if( jt + 1 == jmax){
            /*V velocities on top boundary*/
            V[i][jt] = 0;
            /*U velocities on top boundary (interpolated value)*/
            U[i][jt+1]=2.0-1*U[i][jt];
        }
    }
    
}
