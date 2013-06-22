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
	int t_il, t_ir, t_jb, t_jt, t_omg_i, t_omg_j, t_rank_l, t_rank_r, t_rank_b, t_rank_t;
	int res_i, res_j;
    int quot_i,quot_j;
	MPI_Status status;
    
	MPI_Comm_rank(MPI_COMM_WORLD, myrank) ;
    
	if ( *myrank == 0 ) {
		res_i = ( (imax-1) % iproc) ;
		res_j = ( (jmax-1) % jproc) ;
        quot_i= ( (imax-1) / iproc) ;
        quot_j= ( (jmax-1) / jproc) ;
		*omg_i=1;
		*omg_j=1;
		if( *omg_i <= res_i ){
			*il = (quot_i+1) * ((*omg_i)-1)+1;
			*ir = (quot_i+1) * (*omg_i)+1;
		}
		else{
			*il = ((quot_i+1) * res_i+1)+((*omg_i-res_i-1)*quot_i);
			*ir = ((quot_i+1) * res_i+1)+((*omg_i-res_i)*quot_i);
		}
		if( *omg_j <= res_j ){
			*jb = (quot_j+1) * ((*omg_j)-1)+1;
			*jt = (quot_j+1) * (*omg_j)+1;
		}
		else{
			*jb = ((quot_j+1) * res_j+1)+((*omg_j-res_j-1)*quot_j);
			*jt = ((quot_j+1) * res_j+1)+((*omg_j-res_j)*quot_j);
		}
        
		*rank_l = MPI_PROC_NULL ;
		*rank_r = *myrank + 1 ;
		*rank_t = *myrank + iproc;
		*rank_b = MPI_PROC_NULL ;
        
        
		for (i = 1; i < num_proc; ++i) {
			t_omg_i= (i%iproc) + 1;
			t_omg_j = (i/iproc) + 1;
			

            if( t_omg_i <= res_i ){
                t_il = (quot_i+1) * ((t_omg_i)-1)+1;
                t_ir = (quot_i+1) * (t_omg_i)+1;
            }
            else{
                t_il = ((quot_i+1) * res_i+1)+((t_omg_i-res_i-1)*quot_i);
                t_ir = ((quot_i+1) * res_i+1)+((t_omg_i-res_i)*quot_i);
            }
            if( t_omg_j <= res_j ){
                t_jb = (quot_j+1) * ((t_omg_j)-1)+1;
                t_jt = (quot_j+1) * (t_omg_j)+1;
            }
            else{
                t_jb = ((quot_j+1) * res_j+1)+((t_omg_j-res_j-1)*quot_j);
                t_jt = ((quot_j+1) * res_j+1)+((t_omg_j-res_j)*quot_j);
            }
            
            
            if (t_omg_i == 1){
                t_rank_l = MPI_PROC_NULL ;
                t_rank_r = i + 1 ;
            }
            else if (t_omg_i == iproc ) {
                
                t_rank_r =  MPI_PROC_NULL ;
                t_rank_l = i - 1 ;
            }
            else {
                t_rank_r = i + 1 ;
                t_rank_l = i - 1 ;
			}
            
            if (t_omg_j == 1){
                t_rank_b = MPI_PROC_NULL ;
                t_rank_t = i + iproc;
            }
            else if (t_omg_j == jproc){
                t_rank_t =  MPI_PROC_NULL ;
                t_rank_b = i - iproc;
            }
            else{
                t_rank_t = i + iproc;
                t_rank_b = i - iproc;
			}
            
            
			MPI_Send(&t_omg_i, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&t_omg_j, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
			MPI_Send(&t_il, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
			MPI_Send(&t_ir, 1, MPI_INT, i, 4, MPI_COMM_WORLD);
			MPI_Send(&t_jb, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
			MPI_Send(&t_jt, 1, MPI_INT, i, 6, MPI_COMM_WORLD);
			MPI_Send(&t_rank_l, 1, MPI_INT, i, 7, MPI_COMM_WORLD);
			MPI_Send(&t_rank_r, 1, MPI_INT, i, 8, MPI_COMM_WORLD);
			MPI_Send(&t_rank_b, 1, MPI_INT, i, 9, MPI_COMM_WORLD);
			MPI_Send(&t_rank_t, 1, MPI_INT, i, 10, MPI_COMM_WORLD);
		}
	}
	else{
		MPI_Recv(omg_i, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(omg_j, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
		MPI_Recv(il, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
		MPI_Recv(ir, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);
		MPI_Recv(jb, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
		MPI_Recv(jt, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, &status);
		MPI_Recv(rank_l, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
		MPI_Recv(rank_r, 1, MPI_INT, 0, 8, MPI_COMM_WORLD, &status);
		MPI_Recv(rank_b, 1, MPI_INT, 0, 9, MPI_COMM_WORLD, &status);
		MPI_Recv(rank_t, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
	}
    MPI_Barrier(MPI_COMM_WORLD);
}

void pressure_comm(
                   double **P,
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
    
    int i,j ;
    
    /************************************************/
    /*  send to the left  -  receive form the right */
    /************************************************/
    for (j = jb-1; j <= jt+1; ++j) {
        bufSend[j - jb + 1] = P[il][j] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_l,1,bufRecv,chunk,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD,status) ;
    
    for (j = jb-1; j <= jt+1; ++j) {
        P[ir+1][j] = bufRecv[j-jb + 1] ;
    }
    
    /****************************************/
    /*send to the right  - receive from left*/
    /****************************************/
    
    for (j = jb-1; j <= jt+1; ++j) {
        bufSend[j - jb + 1] = P[ir][j] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_r,2,bufRecv,chunk,MPI_DOUBLE,rank_l,2,MPI_COMM_WORLD,status) ;
    
    for (j = jb-1; j <= jt+1; ++j) {
        P[il-1][j] = bufRecv[j-jb + 1] ;
    }
    /******************************************/
    /* send to the top  - receive from bottom */
    /******************************************/
    
    for (i = il-1; i <= ir+1; ++i) {
        bufSend[i - il + 1] = P[i][jt] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_t,3,bufRecv,chunk,MPI_DOUBLE,rank_b,3,MPI_COMM_WORLD,status) ;
    
    for (i = il-1; i <= ir+1; ++i) {
        P[i][jb-1] = bufRecv[i - il + 1] ;
    }
    
    /**********************************************/
    /* send to the bottom  - receive from the top */
    /**********************************************/
    
    for (i = il-1; i <= ir+1; ++i) {
        bufSend[i - il + 1] = P[i][jb] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_b,4,bufRecv,chunk,MPI_DOUBLE,rank_t,4,MPI_COMM_WORLD,status) ;
    
    for (i = il-1; i <= ir+1; ++i) {
        P[i][jt+1] = bufRecv[i - il + 1] ;
    }
    
    
    
}

void uv_comm(
             double **U,
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
    
    int i,j ;
    
    /************************************************/
    /*  send to the left  -  receive form the right */
    /************************************************/
    
    /* Calculations for U */
    
    for (j = jb-1; j <= jt+1; ++j) {
        bufSend[j - jb + 1] = U[il][j] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_l,1,bufRecv,chunk,MPI_DOUBLE,rank_r,1,MPI_COMM_WORLD,status) ;
    
    for (j = jb-1; j <= jt+1; ++j) {
        U[ir+1][j] = bufRecv[j-jb + 1] ;
    }
    
    /* Calculations for V */
    
    for (j = jb-2; j <= jt+1; ++j) {
        bufSend[j - jb + 2] = V[il][j] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_l,2,bufRecv,chunk,MPI_DOUBLE,rank_r,2,MPI_COMM_WORLD,status) ;
    
    for (j = jb-2; j <= jt+1; ++j) {
        V[ir+1][j] = bufRecv[j-jb + 2] ;
    }
    
    
    
    /****************************************/
    /*send to the right  - receive from left*/
    /****************************************/
    
    /* Calculations for U */
    
    for (j = jb-1; j <= jt+1; ++j) {
        bufSend[j - jb + 1] = U[ir-1][j] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_r,3,bufRecv,chunk,MPI_DOUBLE,rank_l,3,MPI_COMM_WORLD,status) ;
    
    for (j = jb-1; j <= jt+1; ++j) {
        U[il-2][j] = bufRecv[j-jb + 1] ;
    }
    
    /* Calculations for V */
    
    for (j = jb-2; j <= jt+1; ++j) {
        bufSend[j - jb + 2] = V[ir][j] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_r,4,bufRecv,chunk,MPI_DOUBLE,rank_l,4,MPI_COMM_WORLD,status) ;
    
    for (j = jb-2; j <= jt+1; ++j) {
        V[il-1][j] = bufRecv[j-jb + 2] ;
    }
    
    
    /******************************************/
    /* send to the top  - receive from bottom */
    /******************************************/
    
    /* Calculations for U */
    
    for (i = il-2; i <= ir+1; ++i) {
        bufSend[i - il + 2] = U[i][jt] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_t,5,bufRecv,chunk,MPI_DOUBLE,rank_b,5,MPI_COMM_WORLD,status) ;
    
    for (i = il-2; i <= ir+1; ++i) {
        U[i][jb-1] = bufRecv[i - il + 2] ;
    }
    
    /* Calculations for V */
    
    for (i = il-1; i <= ir+1; ++i) {
        bufSend[i - il + 1] = V[i][jt-1] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_t,6,bufRecv,chunk,MPI_DOUBLE,rank_b,6,MPI_COMM_WORLD,status) ;
    
    for (i = il-1; i <= ir+1; ++i) {
        V[i][jb-2] = bufRecv[i - il + 1] ;
    }
    
    
    /**********************************************/
    /* send to the bottom  - receive from the top */
    /**********************************************/
    
    /* Calculations for U */
    
    for (i = il-2; i <= ir+1; ++i) {
        bufSend[i - il + 2] = U[i][jb] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_b,7,bufRecv,chunk,MPI_DOUBLE,rank_t,7,MPI_COMM_WORLD,status) ;
    
    for (i = il-2; i <= ir+1; ++i) {
        U[i][jt+1] = bufRecv[i - il + 2] ;
    }
    
    /* Calculations for V */
    
    for (i = il-1; i <= ir+1; ++i) {
        bufSend[i - il + 1] = V[i][jb] ;
    }
    
    MPI_Sendrecv(bufSend,chunk,MPI_DOUBLE,rank_b,8,bufRecv,chunk,MPI_DOUBLE,rank_t,8,MPI_COMM_WORLD,status) ;
    
    for (i = il-1; i <= ir+1; ++i) {
        V[i][jt+1] = bufRecv[i - il + 1] ;
    }
    
    
}



