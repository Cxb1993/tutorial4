#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void init_parallel(
		int iproc,			/*number of processes in i direction*/
		int jproc,			  /*number of processes in j direction*/
		int imax,			/*maximal number of grid points in i direction*/
		int jmax,			/*maximal number of grid points in j direction*/
		int *myrank,		/*integer denoting ID of the process*/
		int *il,			/*left most point in i direction in subdomain omega */
		int *ir,
		int *jb,
		int *jt,
		int *rank_l,
		int *rank_r,
		int *rank_b,
		int *rank_t,
		int *omg_i,
		int *omg_j,
		int num_proc)
{

	int i ;

MPI_Comm_rank(MPI_COMM_WORLD, &myrank) ;

/*Defining omegas !!!!!!!!!!  */


switch (omg_i) {
	case 1:
		rank_l = MPI_PROC_NULL ;
		rank_r = myrank + 1 ;
		break;
	case iproc:
		rank_r =  MPI_PROC_NULL ;
		rank_l = myrank - 1 ;
		break ;
	default:
		rank_r = myrank + 1 ;
		rank_l = myrank - 1 ;
		break;
}

switch (omg_j) {
	case 1:
		rank_b = MPI_PROC_NULL ;
		rank_t = myrank + iproc - 1  ;
		break;
	case jproc:
		rank_t =  MPI_PROC_NULL ;
		rank_b = myrank  - iproc - 1 ;
		break ;
	default:
		rank_t = myrank - iproc - 1 ;
		rank_b = myrank - iproc - 1 ;
		break;
}




}

void pressure_comm( double **P,
					int il ,
					int ir,
					int jb,
					int jt ,
					int rank_l,
					int rank_r,
					int rank_b,
					int rank_t,
					double *bufSend,
					double *bufRecv,
					MPI_Status *status,
					int chunk )

{


/*put pressure communication here*/



}

void uv_comm(	    double **U,
					double **V,
					int il ,
					int ir,
					int jb,
					int jt ,
					int rank_l,
					int rank_r,
					int rank_b,
					int rank_t,
					double *bufSend,
					double *bufRecv,
					MPI_Status *status,
					int chunk )

{

/*put pressure communication here*/

}



