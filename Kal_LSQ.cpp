#include "mex.h"
//#include "string.h"

// Hassan August 2009
//RLR_c.cpp: no Kalman
//RLR_c2.cpp: with Kalman

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
    if (nrhs!=2)
		mexErrMsgTxt("Wrong number of inputs");
	if (nlhs!=2)
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
    double *dif_I, *dif_I_K;  // dif_I_K is the Kalman Filtered dif_I
    int Dims1[2];
    Dims1[0] = mm; 
    Dims1[1] = nn;
    plhs[0] = mxCreateNumericArray(2,Dims1,mxDOUBLE_CLASS,mxREAL);
	dif_I = (double *)mxGetData(plhs[0]);
    plhs[1] = mxCreateNumericArray(2,Dims1,mxDOUBLE_CLASS,mxREAL);
	dif_I_K = (double *)mxGetData(plhs[1]);
    
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
    //-------------------------------------------------------------------//
    //---------------------- Part 2, Kalman Filter ----------------------//
    double p_k_, p_k; // p_k_ is a priori error, p_k is a posteriori error
    double z_k; // z_k is the measured intensity
    double x_k; // x_k_ & x_k are a-priori and a-posteriori intensities
    double sigma_n_2;
    double mean4var = 0; //mean calculated for variance calculation
    double mean4var2 = 0; //mean^2 calculated for variance calculation
    double s4var; // individual strain samples for var calc
    for (kk = strt_jj+1; kk<end_jj; kk++ )
    {
        s4var = dif_I[11*mm + kk]; // the samples of 11th line
        mean4var += s4var;
        mean4var2 += (s4var*s4var);
    }
    sigma_n_2 = (mean4var2 - mean4var*mean4var/(end_jj-strt_jj))/(end_jj-strt_jj);
    dif_I_v = (double *) malloc( nn * sizeof(double) ); // One dif_k horizantal line 
    for (kk = 0; kk < mm; kk++)
    {
        for (ii = 0; ii < nn; ii++)
            dif_I_v[ii] = dif_I[ii*mm + kk];
        
        dif_I_K[kk] = dif_I_v[0]; // 1st line of the image, no Kalman
        p_k = 0; // initialize
        x_k = dif_I_v[0];
        for (ii = 1; ii < nn; ii++)
        {
            z_k = dif_I_v[ii]; //current point (observation)
            
            p_k_ = p_k + (z_k - dif_I_v[ii-1])*(z_k - dif_I_v[ii-1]); // update a-priori error
            
            x_k = x_k + p_k_ / (p_k_ + sigma_n_2) * (z_k - x_k);
            p_k = sigma_n_2 * p_k_ / (p_k_ + sigma_n_2);
            
            dif_I_K[ii*mm + kk] = x_k;
            
        }
    }
    
    free(dif_I_v);
    
    
}