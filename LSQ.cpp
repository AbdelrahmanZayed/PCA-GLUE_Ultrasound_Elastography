#include "mex.h"
//#include "string.h"

// Hassan August 2009
//RLR_c.cpp: no Kalman
//RLR_c2.cpp: with Kalman

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
    if (nrhs!=2)
		mexErrMsgTxt("Wrong number of inputs");
	if (nlhs!=1)
		mexErrMsgTxt("Wrong number of outputaas");
        
    //---------------------- inputs -----------------------//
    const double *Im; //Im1 is m x n
    int dd;   // the length of the differentiation kernel
       //---------------------- Get inputs --------------------//
	Im = (double *)mxGetData(prhs[0]);
	dd = mxGetScalar(prhs[1]);
    if(dd % 2 == 0)
        mexErrMsgTxt("dd must be odd");
	
	// dimension variables
    int Dims;
    int *Dims_matr;
    
    // Get the information on the first input
	Dims = mxGetNumberOfDimensions(prhs[0]);
	if (Dims!=2)
		mexErrMsgTxt("The number of dimensions of the first parameter must be 2");
	
	Dims_matr = (int *)mxGetDimensions(prhs[0]);
	int mm = Dims_matr[0]; //dimmensions of Im1 and Im2
	int nn = Dims_matr[1];
	if (mm<(dd+1))
		mexErrMsgTxt("Image depth too small compared to differentiation kernel");    
	if (nn<15)
		mexErrMsgTxt("Image width too small");    

    //------------------------------------------------------//
    //---------------------- outputs -----------------------//
    double *dif_I;  // dif_I_K is the Kalman Filtered dif_I
    int Dims1[2];
    Dims1[0] = mm; 
    Dims1[1] = nn;
    plhs[0] = mxCreateNumericArray(2,Dims1,mxDOUBLE_CLASS,mxREAL);
	dif_I = (double *)mxGetData(plhs[0]);
    
    //-------------------------------------------------------------------//
    //---------------------- Part 1, Least Squares ----------------------//
    double *Im_v, *dif_I_v, cov_xy, sum_dd_1, var_x;
    int ii, jj, kk, dd_2, strt_jj, end_jj;
    
    Im_v = (double *) malloc( mm * sizeof(double) ); // One Im vertical line 
    //memset(contD2v_B, 0, lRF * sizeof(double));//backup copy of contD2v for midA to begining of Im
    
    dd_2 = (dd-1)/2;
    strt_jj = dd_2;
    end_jj = mm - dd_2;
    var_x = 3. /dd_2 /(dd_2+1) / (2*dd_2+1); // denuminator
    
    for (ii = 0; ii < nn; ii++)
    {
        for (kk = 0; kk < mm; kk++)
            Im_v[kk] = Im[ii*mm + kk];
        //--------------------------------------------------------//
        //---------------------- 1st sample ----------------------//
        jj = strt_jj; 
        cov_xy = 0;
        sum_dd_1 = 0;
        for (kk = -dd_2; kk < dd_2+1; kk++)
            cov_xy += (kk * Im_v[kk+dd_2]);
        dif_I[ii*mm + jj] = cov_xy * var_x;
        for (kk = 1; kk < dd; kk++)
            sum_dd_1 += Im_v[kk]; //the sum of the current window's numbers (except the 1st number)
        //-----------------------------------------------------------//
        //---------------------- Other samples ----------------------//
        for (jj = strt_jj+1; jj < end_jj; jj++)
        {
            cov_xy += - sum_dd_1 + dd_2 * (Im_v[jj-dd_2-1] + Im_v[jj+dd_2]);
            dif_I[ii*mm + jj] = cov_xy * var_x;
            sum_dd_1 += - Im_v[jj-dd_2] + Im_v[jj+dd_2];
        }
    }
    free(Im_v);
      
    
}