#include "__cf_.h"
#include "mex.h"
void mexFunction ( int nlhs , mxArray * plhs [ ] , int nrhs , const mxArray *
prhs [ ] ) { if ( nrhs != 2 ) mexErrMsgTxt ( "Wrong number of inputs" ) ; if
( nlhs != 1 ) mexErrMsgTxt ( "Wrong number of outputaas" ) ; const double *
Im ; int dd ; Im = ( double * ) mxGetData ( prhs [ 0 ] ) ; dd = mxGetScalar (
prhs [ 1 ] ) ; if ( dd % 2 == 0 ) mexErrMsgTxt ( "dd must be odd" ) ; int
Dims ; int * Dims_matr ; Dims = mxGetNumberOfDimensions ( prhs [ 0 ] ) ; if (
Dims != 2 ) mexErrMsgTxt (
"The number of dimensions of the first parameter must be 2" ) ; Dims_matr = (
int * ) mxGetDimensions ( prhs [ 0 ] ) ; int mm = Dims_matr [ 0 ] ; int nn =
Dims_matr [ 1 ] ; if ( mm < ( dd + 1 ) ) mexErrMsgTxt (
"Image depth too small compared to differentiation kernel" ) ; if ( nn < 15 )
mexErrMsgTxt ( "Image width too small" ) ; double * dif_I ; int Dims1 [ 2 ] ;
Dims1 [ 0 ] = mm ; Dims1 [ 1 ] = nn ; plhs [ 0 ] = mxCreateNumericArray ( 2 ,
Dims1 , mxDOUBLE_CLASS , mxREAL ) ; dif_I = ( double * ) mxGetData ( plhs [ 0
] ) ; double * Im_v , * dif_I_v , cov_xy , sum_dd_1 , var_x ; int ii , jj ,
kk , dd_2 , strt_jj , end_jj ; Im_v = ( double * ) malloc ( mm * sizeof (
double ) ) ; dd_2 = ( dd - 1 ) / 2 ; strt_jj = dd_2 ; end_jj = mm - dd_2 ;
var_x = 3. / dd_2 / ( dd_2 + 1 ) / ( 2 * dd_2 + 1 ) ; for ( ii = 0 ; ii < nn
; ii ++ ) { for ( kk = 0 ; kk < mm ; kk ++ ) Im_v [ kk ] = Im [ ii * mm + kk
] ; jj = strt_jj ; cov_xy = 0 ; sum_dd_1 = 0 ; for ( kk = - dd_2 ; kk < dd_2
+ 1 ; kk ++ ) cov_xy += ( kk * Im_v [ kk + dd_2 ] ) ; dif_I [ ii * mm + jj ]
= cov_xy * var_x ; for ( kk = 1 ; kk < dd ; kk ++ ) sum_dd_1 += Im_v [ kk ] ;
for ( jj = strt_jj + 1 ; jj < end_jj ; jj ++ ) { cov_xy += - sum_dd_1 + dd_2
* ( Im_v [ jj - dd_2 - 1 ] + Im_v [ jj + dd_2 ] ) ; dif_I [ ii * mm + jj ] =
cov_xy * var_x ; sum_dd_1 += - Im_v [ jj - dd_2 ] + Im_v [ jj + dd_2 ] ; } }
free ( Im_v ) ; }
