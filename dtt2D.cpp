/**************************************************************************
 * MEX file to compute 2D discrete trigonometric transforms in double 
 * precision using FFTW. See dtt2D.m for usage notes.
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
    mwSize numelements;    
    const mwSize *dims;
    int NX, NY, numdims;
    double *input_ptr, *output_ptr;
    fftw_plan plan;
    fftw_r2r_kind DTT_type_x, DTT_type_y;
    
    //--------------------------------------------
    // CHECK NUMBER OF INPUTS AND OUTPUTS
    //--------------------------------------------    
    
    //check for proper number of arguments
    if(nrhs!=2) {
        mexErrMsgTxt("Two inputs are required.");
	} else if(nlhs!=1) {
        mexErrMsgTxt("One output is required.");
	}
    
    //--------------------------------------------
    // CHECK AND ALLOCATE DTT TYPE INPUT
    //--------------------------------------------
    
    //check DTT_type input is real and double precision
    if( !(mxIsDouble(prhs[1]) && !mxIsComplex(prhs[1]))) {
        mexErrMsgTxt("Input for DTT_TYPE must be real, and double precision.");
    }
    
    //get pointer to the DTT_TYPE input and check its size
    double * dtt_type_pointer = mxGetPr(prhs[1]);
    mwSize check_el_num = mxGetNumberOfElements(prhs[1]);
    
    //assign DTT types, first cast double input to 32-bit integer for switch, 
    //and then assign the literals
    if (check_el_num == 1){
        
        //assign scalar input to both x and y directions
        switch ( (int) dtt_type_pointer[0] ) {
            case 1: DTT_type_x = FFTW_REDFT00; DTT_type_y = FFTW_REDFT00; break;
            case 2: DTT_type_x = FFTW_REDFT10; DTT_type_y = FFTW_REDFT10; break;
            case 3: DTT_type_x = FFTW_REDFT01; DTT_type_y = FFTW_REDFT01; break;
            case 4: DTT_type_x = FFTW_REDFT11; DTT_type_y = FFTW_REDFT11; break;
			case 5: DTT_type_x = FFTW_RODFT00; DTT_type_y = FFTW_RODFT00; break;
            case 6: DTT_type_x = FFTW_RODFT10; DTT_type_y = FFTW_RODFT10; break;
            case 7: DTT_type_x = FFTW_RODFT01; DTT_type_y = FFTW_RODFT01; break;
            case 8: DTT_type_x = FFTW_RODFT11; DTT_type_y = FFTW_RODFT11; break;   
            default: mexErrMsgTxt("Input for DTT_TYPE must be an integer between 1 and 8.");
        }            
        
    }
    else if (check_el_num == 2){
        
        //assign x-direction
        switch ( (int) dtt_type_pointer[0] ) {
            case 1: DTT_type_x = FFTW_REDFT00; break;
            case 2: DTT_type_x = FFTW_REDFT10; break;
            case 3: DTT_type_x = FFTW_REDFT01; break;
            case 4: DTT_type_x = FFTW_REDFT11; break;
            case 5: DTT_type_x = FFTW_RODFT00; break;
            case 6: DTT_type_x = FFTW_RODFT10; break;
            case 7: DTT_type_x = FFTW_RODFT01; break;
            case 8: DTT_type_x = FFTW_RODFT11; break;            
            default: mexErrMsgTxt("Input for DTT_TYPE must be an integer between 1 and 8.");
        }   
        
        //assign y-direction
        switch ( (int) dtt_type_pointer[1] ) {
            case 1: DTT_type_y = FFTW_REDFT00; break;
            case 2: DTT_type_y = FFTW_REDFT10; break;
            case 3: DTT_type_y = FFTW_REDFT01; break;
            case 4: DTT_type_y = FFTW_REDFT11; break;  
            case 5: DTT_type_y = FFTW_RODFT00; break;
            case 6: DTT_type_y = FFTW_RODFT10; break;
            case 7: DTT_type_y = FFTW_RODFT01; break;
            case 8: DTT_type_y = FFTW_RODFT11; break;              
            default: mexErrMsgTxt("Input for DTT_TYPE must be an integer between 1 and 8.");
        }         
        
    }
    else {
        mexErrMsgTxt("Input for DTT_TYPE must be scalar or length 2.");
    }
    
    //--------------------------------------------
    // CHECK AND ALLOCATE INPUT AND OUTPUT ARRAYS
    //--------------------------------------------
    
    //check the input matrix is real and double precision
    if( !(mxIsDouble(prhs[0]) && !mxIsComplex(prhs[0])) ) {
        mexErrMsgTxt("Input array must be double precision and real.");
    }  

    //check that the input is 2D
    numdims = mxGetNumberOfDimensions(prhs[0]);
    if ( numdims != 2){
        mexErrMsgTxt("Input array must be 2D.");
    }
    
    //get the dimensions of the input array
    dims = mxGetDimensions(prhs[0]);
    numelements = mxGetNumberOfElements(prhs[0]);
    NX = (int) dims[0];
    NY = (int) dims[1];
             
    //create MATLAB output
    output_mat = plhs[0] = mxCreateNumericArray(numdims, dims, mxDOUBLE_CLASS, mxREAL);
    
    //get pointer to input and output arrays
    input_ptr = (double *) mxGetPr(prhs[0]);
    output_ptr = (double *) mxGetPr(output_mat);
    
    //print dimensions of the input array
    //mexPrintf("DTT Type (%d, %d), Array Dimensions %d by %d\n", (int) dtt_type_pointer[0], (int) dtt_type_pointer[0], NX, NY);  
    
    //--------------------------------------------
    // CREATE AND EXECUTE FFTW PLAN
    //--------------------------------------------    
    
    //create FFTW plan
    plan = fftw_plan_r2r_2d(NY, NX, (double *) input_ptr, (double *) output_ptr, DTT_type_y, DTT_type_x, FFTW_ESTIMATE);    
    
    //execute FFTW
    fftw_execute(plan);    
    
    //--------------------------------------------
    // CLEANUP ALLOCATED VARIABLES
    //--------------------------------------------

    //cleanup
    fftw_destroy_plan(plan);
    
    return;
}