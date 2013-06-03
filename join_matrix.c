#include "helper.h"
#include "limits.h"

/* -------------------------------------------------------------------------------- 

  The program allows one to assemble partial matrices corresponding to different
  parts of a domain to an individual matrix for the whole domain. The matrices 
  to be assembled have to be stored using write_matrix from helper.c.

  The dimension of the resulting matrix is determined by the dimensions of partial
  matrices. The overlapping positions are taken into account.

  Note that x- and y-dimensions of the partial matrices must concide.

  -------------------------------------------------------------------------------- */

#define FILE_ERROR( szFilename, szError ) \
	{ \
	   char szBuff[80]; \
	   sprintf( szBuff, szError, szFilename ); \
	   ERROR( szBuff ); \
	}


int main( int argc, char *argv[] ) 
{

    
  
    int i,j ; 
    double **ResultingMatrix = NULL;

    int imin = INT_MAX;
    int jmin = INT_MAX;
    int imax = INT_MIN;
    int jmax = INT_MIN;
    int nFiles = argc-2;	/* number of matrices to be joined */

    double xlength;
    double ylength;
    double ftmp;


    FILE **ArrFiles;		/* array with file handles*/
    int  *ArrImin;
    int  *ArrImax;
    int  *ArrJmin;
    int  *ArrJmax;

    char str[100];              /* string for fgets  */
    
    int nCount = 0;		/* number of iterations */

    if( nFiles < 1 )
    {
	printf("usage: %s <Matrix_1> <Matrix_2> ... <Matrix_n> <Ausgabefile>\n", argv[0] );
	exit(1);
    }  

    ArrFiles     = (FILE **)malloc( nFiles * sizeof( FILE *));
    ArrImin = (int *)malloc( nFiles * sizeof( int ));
    ArrImax = (int *)malloc( nFiles * sizeof( int ));
    ArrJmin = (int *)malloc( nFiles * sizeof( int ));
    ArrJmax = (int *)malloc( nFiles * sizeof( int ));

    printf("-------------------------------------------------------------------------------\n");
    /* Iterate over all files and generate file handles and */
    /* dimensions of the matrices */
    for( i = 0; i < nFiles; i++)
    {
	const char *szFilename = argv[i+1];		  /* argv[0] is the programm name, argv[1] is the first parameter  */
	FILE *fh = fopen( szFilename, "r" );
	
	if( fh == NULL)  
	    FILE_ERROR( szFilename, "Error: can't open the file %s" ); 
	
	/* xlength */
	if( fscanf( fh, "%lf\n", &ftmp) == 0 )  	    
	    FILE_ERROR( szFilename, "Error: can't read the file %s" ); 
	
	if( i == 0 )  xlength = ftmp;			  /* store xlength of the first file */
	if( xlength != ftmp )
	    FILE_ERROR( szFilename, "Error: wrong value of xlength in the input file %s" ); 
	
	/* ylength */
	if( fscanf( fh, "%lf\n", &ftmp) == 0 )  
	    FILE_ERROR( szFilename, "Error: can't read the file %s" ); 
	
	if( i == 0 )  ylength = ftmp;			  /* store ylength of the first file */
	if( ylength != ftmp )
	    FILE_ERROR( szFilename, "Error: wrong value of ylength in the input file %s" ); 

	/* imin */
	if( fscanf( fh, "%d\n", &(ArrImin[i])) == 0 )  
	    FILE_ERROR( szFilename, "Error: can't read the file %s"); 

	/* imax */
	if( fscanf( fh, "%d\n", &(ArrImax[i])) == 0 )  
	    FILE_ERROR( szFilename, "Error: can't read the file %s" ); 

	/* jmin */
	if( fscanf( fh, "%d\n", &(ArrJmin[i])) == 0 )  
	    FILE_ERROR( szFilename, "Error: can't read the file %s" ); 

	/* jmax */
	/* fscanf ignores new line character. this leads to problems with fread if 
           the first byte is \n.
        /* Therefore we use fgets to read the file linewise */
        if( fgets(str,100,fh) == 0 )
            FILE_ERROR( szFilename, "Error: can't read the file %s" ); 
        if( sscanf( str, "%d\n", &(ArrJmax[i])) == 0 )  
            FILE_ERROR( szFilename, "Error: can't read the file %s" );

	/* calculate the maximal dimension of matrixes */
	if( ArrImin[i] < imin )  imin = ArrImin[i];
	if( ArrImax[i] > imax )  imax = ArrImax[i];
	if( ArrJmin[i] < jmin )  jmin = ArrJmin[i];
	if( ArrJmax[i] > jmax )  jmax = ArrJmax[i];

	printf("name: %s    size:  %5d %5d  %5d %5d \n",
	       szFilename,
	       ArrImin[i],
	       ArrImax[i],
	       ArrJmin[i],
	       ArrJmax[i] );
	
	ArrFiles[i] = fh;				  /* store the file handle in the array */
    }    

    printf("-------------------------------------------------------------------------------\n");

    /* generate the matrix of an appropriate dimension */

    
    ResultingMatrix = matrix( imin, imax, jmin, jmax );

    /* read and write until the end of the first file is reached  */
    while(1) //(j=0;j<nFiles;j++)
    {
	/* initialize matrix  */
	init_matrix( ResultingMatrix, imin, imax, jmin, jmax, 0.0 );   /* TODO NAN */

	/*  read matrix  */
	for( i = 0; i < nFiles; i++)
	{
	    FILE *fh = ArrFiles[i];
	   
           
	    if( read_matrix( fh,
			     ResultingMatrix, 
			     ArrImin[i],
			     ArrImax[i],
			     ArrJmin[i],
			     ArrJmax[i] ) == EOF) 
	    {
		/* no input - no result, exit the programm  */
	          
	      
	     return 0; 
	    }
	}

        
	/* write the assembled matrix without range information*/
	write_matrix( argv[nFiles + 1], 
		      ResultingMatrix, 
		      imin, 
		      imax, 
		      jmin, 
		      jmax, 
		      xlength, 
		      ylength,
		      nCount == 0,1 );

	++nCount;
	printf("%d. iteration\n", nCount);
    }



    return 1;			/* Dummy */
}                

