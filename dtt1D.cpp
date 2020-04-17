/**************************************************************************
 * MEX file to compute 1D discrete trigonometric transforms in double 
 * precision using FFTW. See dtt1D.m for usage notes.
 *
 * author: Bradley Treeby
 * date: 31 May 2012
 * last update: 19 March 2020
 *
 * Copyright (C) 2012-2020 Bradley Treeby
 *
 * This program is free software: you can redistribute it and/or modify it 
 * under the terms of the GNU General Public License as published by the 
 * Free Software Foundation, either version 3 of the License, or (at your 
 * option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program.  If not, see <https://www.gnu.org/licenses/>.
 *************************************************************************/

#include <matrix.h>
#include <mex.h>
#include "fftw3.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    //--------------------------------------------
    // DECLARE VARIABLES
    //--------------------------------------------
    
    mxArray *output_mat;
    const mwSize *dims;
    mwSize check_el_num;
    mwSize numelements;
    int NX, NY, N, numdims;
    int i;
    double *input_ptr, *output_ptr;
    fftw_plan plan;
    int DIM = 1;    // set default DIM to 1 if not given by user
    int rank = 1;   // DTTs are all 1D or rank 1
    int n[1];       // DTT size is thus also rank 1
    int howmany;
    int idist, odist;
    int istride, ostride;
    
    //--------------------------------------------
    // CHECK NUMBER OF INPUTS AND OUTPUTS
    //--------------------------------------------
    
    //check for proper number of input and output arguments
    if( !(nrhs == 2 || nrhs == 3) ) {
        mexErrMsgTxt("Two or three inputs are required.");
	} else if(nlhs!=1) {
        mexErrMsgTxt("One output is required.");
	}
    
    //--------------------------------------------
    // CHECK AND ALLOCATE DTT TYPE INPUT
    //--------------------------------------------
    
    //check DTT_type input is real, scalar, and double precision
    check_el_num = mxGetNumberOfElements(prhs[1]);
    if( !(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) && check_el_num == 1) ){
        mexErrMsgTxt("Input for DTT_TYPE must be real, scalar, and double precision.");
    }
    
    //get the DTT type (cast double input to 32-bit integer)
    int DTT_type = (* (double *) mxGetPr(prhs[1]));
    
    //check the value
    if ( !((DTT_type == 1) || (DTT_type == 2) || (DTT_type == 3) || (DTT_type == 4) || (DTT_type == 5) || (DTT_type == 6) || (DTT_type == 7) || (DTT_type == 8)) ){
        mexErrMsgTxt("Input for DTT_TYPE must be an integer between 1 and 8.");
    }
    
    //--------------------------------------------
    // CHECK AND ALLOCATE DIM INPUT
    //--------------------------------------------
    
    //if DIM input is given, check it
    if (nrhs == 3){
        
        //get the number of elements
        check_el_num = mxGetNumberOfElements(prhs[2]);
        
        //check DIM input is real, scalar, and double precision
        if( !(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]) && check_el_num == 1) ){
            mexErrMsgTxt("Input for DIM must be real, scalar, and double precision.");
        }
        
        //get the value of the DIM input (cast double input to 32-bit integer)
        // (* something) means take the first element == something[0]
        DIM = (int)(* (double *) mxGetPr(prhs[2]));
        
        //check the value is 1 or 2
        if ( !( (DIM == 1) || (DIM == 2) ) ){
            mexErrMsgTxt("Input for DIM must be 1 or 2.");
        }
        
    }
    
    //--------------------------------------------
    // CHECK AND ALLOCATE INPUT AND OUTPUT ARRAYS
    //--------------------------------------------
    
    //check the input matrix is real and double precision
    if( !(mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0])) ) {
        mexErrMsgTxt("Input array must be double precision and real.");
    }  

    //check that the input is either 1D or 2D
    numdims = mxGetNumberOfDimensions(prhs[0]);     //number of dimensions is always >= 2
    if ( numdims > 2){
        mexErrMsgTxt("Input array must be 1D or 2D.");
    }
    
    //get the dimensions of the input array
    dims = mxGetDimensions(prhs[0]);
    numelements = mxGetNumberOfElements(prhs[0]);
    NX = (int)dims[0];
    NY = (int)dims[1];
    N = NX*NY;
    
    //check if 1D and force the correct DIM (input DIM is not used)
    if (NX == 1) {
        DIM = 2;
    }
    else if (NY == 1) {
        DIM = 1;
    }
    
    //print dimensions of the input array
    //mexPrintf("DTT Type %d, Array Dimensions %d by %d, DTT on Dimension %d\n", DTT_type, NX, NY, DIM);
   
    //create MATLAB output
    output_mat = plhs[0] = mxCreateNumericArray(numdims, dims, mxDOUBLE_CLASS, mxREAL);
    
    //get pointer to input and output arrays
    input_ptr = (double *) mxGetPr(prhs[0]);
    output_ptr = (double *) mxGetPr(output_mat);
    
    //--------------------------------------------
    // DEFINE PLAN VARIABLES
    //--------------------------------------------
    
    //define plan_many variables based on DIM setting (these are switched 
    //to c++ as MATLAB inputs are column major)
    switch ( DIM ) {
        case 1 :
            //perform DTT over columns of input array
            n[0] = NX;
            howmany = NY;
            idist = odist = NX;
            istride = ostride = 1;            
            break;
        case 2 :      
            //perform DTT over rows of input array
            n[0] = NY;
            howmany = NX;
            idist = odist = 1;
            istride = ostride = NX;
            break;
        default :
            mexErrMsgTxt("Input for DIM must be 1 or 2.");
    }
    
    //plan_many requires a pointer to the input type. However, FFTW_REDFT00 etc
    //are defined using #define so we cannot directly get a pointer to them. 
    //Instead we allocate some memory, copy the value of the input type, and use the 
    //pointer to this memory.
    fftw_r2r_kind * fKind = (fftw_r2r_kind *) fftw_malloc(sizeof(fftw_r2r_kind));
     
    //set DTT type
    switch ( DTT_type ) {
        case 1: fKind[0] = FFTW_REDFT00; break;
        case 2: fKind[0] = FFTW_REDFT10; break;
        case 3: fKind[0] = FFTW_REDFT01; break;
        case 4: fKind[0] = FFTW_REDFT11; break;
        case 5: fKind[0] = FFTW_RODFT00; break;
        case 6: fKind[0] = FFTW_RODFT10; break;
        case 7: fKind[0] = FFTW_RODFT01; break;
        case 8: fKind[0] = FFTW_RODFT11; break;  		
        default: mexErrMsgTxt("Input for DTT_TYPE must be an integer between 1 and 8.");
    }
    
    //--------------------------------------------
    // CREATE AND EXECUTE FFTW PLAN
    //--------------------------------------------
    
    /*fftw_plan fftw_plan_many_r2r(int rank, const int *n, int howmany,
                                  double *in, const int *inembed,
                                  int istride, int idist,
                                  double *out, const int *onembed,
                                  int ostride, int odist,
                                  const fftw_r2r_kind *kind, unsigned flags);*/
    
    //create plan using the input and output pointers directly (out of place transform)
    plan = fftw_plan_many_r2r(rank, (int *) &n, howmany, (double *) input_ptr, NULL, istride, idist, (double *) output_ptr, NULL, ostride, odist, (fftw_r2r_kind*) fKind, FFTW_ESTIMATE);    

    //execute plan
    fftw_execute(plan);    
    
    //--------------------------------------------
    // CLEANUP ALLOCATED VARIABLES
    //--------------------------------------------
    
    //cleanup
    fftw_destroy_plan(plan);
    fftw_free((fftw_r2r_kind*)fKind);
    
    return;
}