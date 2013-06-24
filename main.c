#include "helper.h"
#include "stddef.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include "parallel.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/* CFD Lab - Worksheet 4 - Group 3
 * Camacho Barranco, Roberto
 * Gavranovic, Stefan
 * Valizadeh, Mahyar
 */

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
int main(int argc, char** argv){
    /*Geometry data*/
    double xlength;		 /*domain size in x-direction*/
    double ylength;		/*domain size in y-direction*/
    /* The computational domain is thus = [0; xlaenge] X [0; ylaenge].*/
    int imax;		/* number of interior cells in x-direction*/
    int jmax;		/* number of interior cells in y-direction*/
    double dx;		/* length dx of one cell in x-direction*/
    double dy;		/* length dy of one cell in y-direction*/
    /* Time-stepping data:*/
    double t;		/* current time value*/
    double t_end;		/* final time tend*/
    double dt;		/* time step size dt*/
    double tau;		/* safety factor for time step size control T*/
    double dt_value;	/* time interval for writing visualization data in a file*/
    int n;			/* current time iteration step*/
    /* Pressure iteration data:*/
    int itermax;		/* maximum number of pressure iterations in one time step*/
    int it;			/* SOR iteration counter*/
    double res;		/* residual norm of the pressure equation*/
    double eps;		/* accuracy criterion epsilon (tolerance) for pressure iteration (res < eps)*/
    double omg;		/* relaxation factor omega for SOR iteration*/
    double alpha;		/* upwind differencing factor alpha (see equation (4))*/
    /* Problem-dependent quantities:*/
    double Re;		/* Reynolds number Re*/
    double GX,GY;		/* external forces gx; gy, e.g. gravity*/
    double UI,VI,PI;	/* initial data for velocities and pressure*/
    /* Arrays*/
    double **U;		/* velocity in x-direction*/
    double **V;		/* velocity in y-direction*/
    double **P;		/* pressure*/
    double **RS;		/* right-hand side for pressure iteration*/
    double **F,**G;     /* F;G*/
    double *bufSend;
    double *bufRecv ;
    int n_div;
    int iproc;
    int jproc;
    int myrank;
    int il;
    int ir;
    int jb;
    int jt;
    int rank_l;
    int rank_r;
    int rank_b;
    int rank_t;
    int omg_i;
    int omg_j;
    int num_proc;
    int chunk;
    MPI_Status *status;
    int a,b ;
    int i,j;
    
    MPI_Init( &argc, &argv );                    /* execute n processes      */
    MPI_Comm_size( MPI_COMM_WORLD, &num_proc);     /* asking for the number of processes  */
    
    /* read the program configuration file using read_parameters()*/
    read_parameters("cavity100.dat", &Re, &UI, &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau, &itermax, &eps, &dt_value, &iproc, &jproc);
    
    init_parallel(iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t, &omg_i, &omg_j, num_proc);
    
    /*calculating max for chunk size*/
    
    a = jb - jt + 4 ;
    b = ir - il + 4 ;
    
    if (a>b) chunk = a;
    else chunk = b ;
    
    bufSend = (double *)  malloc((size_t)( chunk * sizeof( double )));
    bufRecv = (double *)  malloc((size_t)( chunk * sizeof( double )));
    
    /* set up the matrices (arrays) needed using the matrix() command*/
    U = matrix(il-2, ir+1, jb-1, jt+1);
    V = matrix(il-1, ir+1, jb-2, jt+1);
    P = matrix(il-1, ir+1, jb-1, jt+1);
    RS = matrix(il, ir, jb, jt);
    F = matrix(il-2, ir+1, jb-1, jt+1);
    G = matrix(il-1, ir+1, jb-2, jt+1);
    
    
    /* initialize current time and time step*/
    t = 0;
    n = 0;
    
    /* create the initial setup init_uvp()*/
    init_uvp(UI, VI, PI,  il, ir, jb, jt, U, V, P);
    
    if (myrank==0){
    write_matrix("matrix.dat",U,il,ir+1,jb,jt+1,xlength,ylength,1,0);
    write_matrix("matrix.dat",V,il,ir+1,jb,jt+1,xlength,ylength,0,0);
    }
    
    
    /* ----------------------------------------------------------------------- */
    /*                             Performing the main loop                    */
    /* ----------------------------------------------------------------------- */
    
    
    while (t<=t_end){
        
        calculate_dt(Re, tau, &dt, dx, dy, il, ir, jb, jt, U, V, num_proc, myrank);
        boundaryvalues(il, ir, jb, jt, imax, jmax, U, V);
        calculate_fg(Re, GX, GY, alpha, dt, dx, dy, il, ir, jb, jt, imax, jmax, U, V ,F , G);
        calculate_rs(dt, dx, dy, il, ir, jb, jt, F, G, RS);
        res = 1.0;
        it = 0;
        
        while(it < itermax && res > eps){
            sor( omg, dx, dy, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, bufSend, bufRecv, status, chunk , P, RS, myrank, imax, jmax,&res,omg_i,omg_j,iproc,jproc);
            if (myrank==0){
/*                printf("residual=%f\n",res);*/
            }
            it++;
        }
        
        
        
        calculate_uv(dt, dx, dy, il, ir, jb, jt, U, V, F, G, P);
        
        uv_comm(U,V,il,ir,jb,jt,rank_l,rank_r,rank_b,rank_t,bufSend, bufRecv, status, chunk);
        
        n_div=(dt_value/dt);
        if (n % n_div == 0) {
            output_uvp( U, V, P, il, ir, jb, jt, omg_i, omg_j,"cavity",n,dx,dy, myrank,imax,jmax);
        }
        
        t = t + dt;
        n++;
    }
    
    
    
    if (myrank==0){
        write_matrix("matrix.dat",U,il,ir+1,jb,jt+1,xlength,ylength,0,0);
        
        write_matrix("matrix.dat",V,il,ir+1,jb,jt+1,xlength,ylength,0,0);
    }
    
    
    
    free(bufSend);
    free(bufRecv);
    /* Destroy memory allocated*/
    free_matrix(U, il-2, ir+1, jb-1, jt+1);
    free_matrix(V, il-1, ir+1, jb-2, jt+1);
    free_matrix(P, il-1, ir+1, jb-1, jt+1);
    free_matrix(RS, il, ir, jb, jt);
    free_matrix(F, il-2, ir+1, jb-1, jt+1);
    free_matrix(G, il-1, ir+1, jb-2, jt+1);
    Programm_Stop("finished its work");
    
    return 0;
}
