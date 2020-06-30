/*
 * File condlike.cpp
 *
 * Part of the `epiPOMS' R package
 *
 * Calculation of the conditional likelihood required by iFFBS
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    
    SEXP condlike(SEXP len1, SEXP R1, SEXP F1, SEXP X1, SEXP Types1, SEXP thetaR1, SEXP thetaF1){
        int *R, *F, *X;
        double *thetaR, *thetaF;
        int  rf, x, j1, len, Types;
        SEXP N, Nrx, Nfx;
        
        R1 = coerceVector(R1, INTSXP);
        F1 = coerceVector(F1, INTSXP);
        X1 = coerceVector(X1, INTSXP);
        len1 = coerceVector(len1, INTSXP);
        Types1 = coerceVector(Types1, INTSXP);
        thetaR1 = coerceVector(thetaR1, REALSXP);
        thetaF1 = coerceVector(thetaF1, REALSXP);
        
        R = INTEGER(R1);
        F = INTEGER(F1);
        X = INTEGER(X1);
        Types = INTEGER(Types1)[0];
        len = INTEGER(len1)[0];
        thetaR = REAL(thetaR1);
        thetaF = REAL(thetaF1);
        
        PROTECT(Nrx = allocVector(REALSXP, len) );
        PROTECT(Nfx = allocVector(REALSXP, len) );
        PROTECT(N = allocVector(REALSXP, len ));
        
        //  set to zero
        for (j1=0 ; j1<len ; j1++) {
            REAL(Nrx)[j1] = 0.0;
            REAL(Nfx)[j1] = 0.0;
        }
        
        for (j1=0; j1<len; j1++){
            for (x=0; x<Types; x++ ){
                for (rf = 0; rf<(Types+1); rf++){
                    if ( (X[j1]==x) && (R[j1]==rf)) {
                        REAL(Nrx)[j1] = thetaR[x*(Types+1) + rf];
                    }  
                    if ( (X[j1]==x) && (F[j1]==rf)) {
                        REAL(Nfx)[j1] = thetaF[x*(Types+1) + rf];
                    }  
                }   
            }
        }
        
        for (j1=0; j1<len; j1++){
            REAL(N)[j1] = REAL(Nrx)[j1] * REAL(Nfx)[j1];
        }
        
        UNPROTECT(3);
        return N;
    }        
}      
