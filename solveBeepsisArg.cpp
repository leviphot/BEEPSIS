/*
    This file is part of the BEEPSIS toolbox.
    See LICENSE.md for information about using/distributing this file.

    Functions in this files are to be called internally by the beepsis.m

    Compile by the mex compiler
      mex -R2018a solveBeepsisArg.cpp

*/

#include <cmath>
#include <iostream>
#include <cstring>

#include "mex.h"

using namespace std;


void symbandchol(const int N, const int B, double *C)
{
    for (int kk = 0; kk<N; ++kk)
    {
        int last = (kk+B-1<N) ? B : N - kk + 1;
        for (int jj = 1; jj<last; ++jj)
        {
            int ii = kk + jj;
            double ratio = C[jj+kk*B] / C[kk*B];
            for (int mm = 0; mm<last-jj; ++mm)
                C[mm + ii*B] -= ratio* C[mm+jj + kk*B];            
        }
        double sq = sqrt(C[kk*B]);
        for (int jj = 0; jj<B; ++jj)
            C[jj + kk*B] /= sq;
    }
    for (int is = (N-B+2)*B-1; is<N*B; is+=B-1)
        for (int ll = is; ll<N*B; ll+=B)
            C[ll] = 0;    
}

double solve_sb_rhs(const int N, const int B, const double *CH, const double * rs, double *x = NULL)
{
    bool isxlocal = x == NULL;
    if (isxlocal)
        x = new double[N];

    x[N-1] = rs[N-1] / CH[(N-1)*B];
    for (int kk = N-2; kk>=0 ; kk--)
    {
        int last = (N<kk+B) ? N : kk+B;
        double sm = 0;
        for (int ll = kk+1, ic=1; ll<last; ++ll, ++ic)
            sm += CH[ic + kk*B] * x[ll];
        x[kk] = (rs[kk]-sm)/CH[kk*B];
    }
    
    x[0] /= CH[0];
    for(int kk = 1; kk<N; kk++)
    {
        int last = B<kk+1 ? B : kk+1;
        double sm = 0;
        for (int ll = 1; ll<last; ll++)
            sm += CH[ll + (kk-ll)*B] * x[kk-ll];
        x[kk] -= sm;
        x[kk] /= CH[kk*B];
    }
    double r = 0;
    for (int kk = 0; kk<N; kk++)
        r += x[kk]*rs[kk];
    
    if (isxlocal)
        delete [] x;
    
    return r;
}

double logdetC(const int N, const int B, const double *CH)
{
    double r = 0;
    for (int kk = 0; kk<N; kk ++)
        r += log(CH[kk*B]);
    return 2*r;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
/* Argument checks */
    if( nlhs > 4 ) {
        mexErrMsgTxt("Too many outputs");
    }
    if( nrhs != 2 ) {
        mexErrMsgTxt("Need exactly 2 inputs");
    }
    
    //cov matrix
    if( !(mxIsDouble(prhs[0]) ) || mxIsSparse(prhs[0]) || mxIsComplex(prhs[0]) || (mxGetNumberOfDimensions(prhs[0]) > 2))
    {
        mexErrMsgTxt("1st argument must be real double full matrix");
    }
    //rhs vector/matrix
    if( !(mxIsDouble(prhs[1]) ) || mxIsSparse(prhs[1]) || mxIsComplex(prhs[1]) || (mxGetNumberOfDimensions(prhs[1]) > 2) )
    {
        mexErrMsgTxt("2nd argument must be real double full vector/matrix");
    }
    
    const mwSize *sz = mxGetDimensions(prhs[0]);
    const int B = sz[0];
    const int N = sz[1];
    
    if (N<1)
        mexErrMsgTxt("Symmetric band array size must be at least 1");
    if (B<1)
        mexErrMsgTxt("Symmetric band array number of bands must be at least 1");
    if (B>N)
        mexErrMsgTxt("Symmetric band array number of bands must <= array size");


    const mwSize *sz2 = mxGetDimensions(prhs[1]); 
    if (sz2[0] != N)
    {
        mexErrMsgTxt("Symmetric band array size must be the same as number of rows of rhs");
    }
    
    const int RHC = sz2[1];
    
        
    mwSize ss[2]; 
    ss[0] = 1; ss[1] = 1;
    mxArray * mLogDet = mxCreateNumericArray(2, ss, mxDOUBLE_CLASS, mxREAL);
    
    ss[1] = RHC;
    mxArray * mSol    = mxCreateNumericArray(2, ss, mxDOUBLE_CLASS, mxREAL);
    
    ss[0] = B; ss[1] = N;
    mxArray *mChol   = mxCreateNumericArray(2, ss, mxDOUBLE_CLASS, mxREAL);
    
    ss[0] = N; ss[1] = 1;
    mxArray *mVec    = (nlhs<4) ? NULL: mxCreateNumericArray(2, ss, mxDOUBLE_CLASS, mxREAL);;
    
    
    const double *bandarray = mxGetDoubles(prhs[0]);
    const double *rharray  = mxGetDoubles(prhs[1]);
    double *sol    = mxGetDoubles(mSol);
    double *logdet = mxGetDoubles(mLogDet);
    double *chol   = mxGetDoubles(mChol);
    memcpy(chol, bandarray, N*B*sizeof(double));
        
    //calculate Choleski decomposition of covariance matrix
    symbandchol(N, B, chol);
    //calculate log det of covariance matrix
    logdet[0] = logdetC(N, B, chol);
    
    //calculate 1st sum of mC^{-1}m'
    sol[0]    = solve_sb_rhs(N, B, chol, &rharray[0], mVec ? mxGetDoubles(mVec) : NULL);
    
    //calculate other sums mC^{-1}m'
    for (int kk = 1; kk<RHC; kk++)
        sol[kk]    = solve_sb_rhs(N, B, chol, &rharray[kk*N]);

    
    // return results or destroy unused variables
    if (nlhs>=1)
        plhs[0] = mSol;
    else
        mxDestroyArray(mSol);
    
    if (nlhs >= 2)
        plhs[1] = mLogDet;
    else
        mxDestroyArray(mLogDet);
    
    if (nlhs >= 3)
        plhs[2] = mChol;
    else
        mxDestroyArray(mChol);
    
    if (nlhs >= 4)
        plhs[3] = mVec;
    
}