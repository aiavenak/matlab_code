/*
 *  box_count_3d.cpp
 *  
 *
 *  Created by andrujalabruja on 9/26/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <cmath>    //  Definitions for Matlab API
#include <mex.h>    //  Definitions for Matlab API

using namespace std;

void box_count_3d(double o[],double N[], int r, int dim[])
{
	int d1 = dim[0];
	int d2 = dim[1];
	int d3 = dim[2];
	int cont = 0, n = 0;
	int i, j, k;
	int ii, jj, kk;
	
	// printf("d1,d2,d3: %d, %d, %d \n",d1,d2,d3);
	// printf("r: %d\n",r);
	
	// loop over N
	for (i = 0;i < d1/r; i++) {
		for (j = 0;j < d2/r; j++) {
			for (k = 0;k < d3/r; k++) {
				
				//loop over o
				for (ii = i*r;ii < (i+1)*r; ii++) {
					for (jj = j*r;jj < (j+1)*r; jj++) {
						for (kk = k*r;kk < (k+1)*r; kk++) {
							n = n + o[ii + jj*d1 + kk*d1*d2]; 
						}
					}
				}	
				//return result in pN[]
				// printf("n: %d\n",n);
				N[i + j*d1/r + k*(d1/r)*(d2/r)] = n;
				// printf("N[%d]: %d\n",i + j*d1/r + k*(d1/r)*(d2/r),N[i + j*d1/r + k*(d1/r)*(d2/r)]);
				n = 0;
				
			}
		}
	}
	
}

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[] )
{
	
	double *pr, *po, *pN, *pdim_vect;
	int i , r;
	int dim_vect_N[3], dim_vect_o[3];
	const int *pdim;
	
	/* Check for proper number of input and output arguments. */
	if (nrhs != 2) {
		mexErrMsgTxt("Two input arguments required.");
	}
	if (nlhs > 1) {
		mexErrMsgTxt("Only one output argument is allowed.");
	}
	
	/* get numeber of dimensions and the size of each dimension */
	//num_dim = mxGetNumberOfDimensions(prhs[0]);
	pdim = mxGetDimensions(prhs[0]);
	
	/* size of sliding windlow*/
	
	pr = mxGetPr(prhs[1]);
	r = (int)(*pr);
	//r = (int)(pr[0]);
	
	/* Allocate the space for the return argument.*/
	/* dimensions of output are dim_vect_o[i]/r */
	for (i = 0;i<3;i++) {
		dim_vect_o[i] = (int)(*pdim++);
		dim_vect_N[i] = (int)(dim_vect_o[i]/(r));
	}
	
	plhs[0] = mxCreateNumericArray(3, dim_vect_N,mxDOUBLE_CLASS, mxREAL);
	
	
        // Assign pointers to each input and output
        po = mxGetPr(prhs[0]);
	pN = mxGetPr(plhs[0]);
	
	// Call the subroutine
	box_count_3d(po, pN, r, dim_vect_o);
}










