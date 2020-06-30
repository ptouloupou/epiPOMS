/*
 * File statetrans.cpp
 *
 * Part of the `epiPOMS' R package
 *
 * Calculation of total number of number of state-specific transitions over the entire time 
 * period
 */


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    SEXP statetrans(SEXP t1, SEXP i1, SEXP X1){
        int *X;
        int  T, I, i, t, p;
        SEXP N;
        
        
        X1 = coerceVector(X1, INTSXP);
        t1 = coerceVector(t1, INTSXP);
        i1 = coerceVector(i1, INTSXP);
        
        X = INTEGER(X1);
        T = INTEGER(t1)[0];
        I = INTEGER(i1)[0];
        
        PROTECT(N = allocVector(INTSXP, 5) );
        
        //  set to zero
        for (p=0; p<5; p++){
            INTEGER(N)[p] = 0;
        }
        
        
        for (t=1; t<T; t++ ){
            for (i =0 ; i<I ; i++) {
                if ((X[(t-1)*I+i]>0) && (X[t*I+i]>0) && (X[(t-1)*I+i]!= X[t*I+i]) ) {
                    INTEGER(N)[0]  += 1;
                }
                if ((X[(t-1)*I+i]==0) && (X[t*I+i]>0) ) {
                    INTEGER(N)[1]  += 1;
                }
                if ((X[(t-1)*I+i]>0) && (X[t*I+i]==0) ) {
                    INTEGER(N)[2]  += 1;
                }
                if ((X[(t-1)*I+i]==0) && (X[t*I+i]==0) ) {
                    INTEGER(N)[3]  += 1;
                }
                if ((X[(t-1)*I+i]>0) && (X[t*I+i]>0) && (X[(t-1)*I+i]== X[t*I+i]) ) {
                    INTEGER(N)[4]  += 1;
                }
            }
        }
        
        
        UNPROTECT(1);
        return N;
    }
}
