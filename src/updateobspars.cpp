/*
 * File updateobspars.cpp
 *
 * Part of the `epiPOMS' R package
 *
 * Calculations required to update the observation parameters
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    SEXP updateobspars(SEXP len1, SEXP R1, SEXP F1, SEXP X1, SEXP Types1){
        int *R, *F, *X;
        int  j, len, Types;
        SEXP Counts;
        
        R1 = coerceVector(R1, INTSXP);
        F1 = coerceVector(F1, INTSXP);
        X1 = coerceVector(X1, INTSXP);
        len1 = coerceVector(len1, INTSXP);
        Types1 = coerceVector(Types1, INTSXP);
        
        R = INTEGER(R1);
        F = INTEGER(F1);
        X = INTEGER(X1);
        Types = INTEGER(Types1)[0];
        len = INTEGER(len1)[0];
        
        
        PROTECT(Counts = allocVector(INTSXP, 9));
        
        //  set to zero
        for (j=0 ; j<9 ; j++) {
            INTEGER(Counts)[j] = 0;
        }
        
        for (j=0; j<len; j++){
            
            if ( X[j] > 0 ) {
                
                if ( R[j] == 0 ) { INTEGER(Counts)[1] +=1;
                } else {
                    INTEGER(Counts)[0] +=1;
                    if ( R[j] < Types + 1 ){
                        if ( R[j] == Types ) {
                            if (X[j] == Types) { INTEGER(Counts)[8] +=1;
                            } else {INTEGER(Counts)[6] +=1; }
                        } else {
                            if (X[j] == Types) { INTEGER(Counts)[7] +=1;
                            } else {
                                if (X[j] == R[j]) { INTEGER(Counts)[4] +=1;
                                } else {INTEGER(Counts)[5] +=1; }
                                
                            }
                        }
                    }
                }
                
                if ( F[j] == 0 ) { INTEGER(Counts)[3] +=1;
                } else {
                    INTEGER(Counts)[2] +=1;
                    if ( F[j] < Types + 1 ){
                        if ( F[j] == Types ) {
                            if (X[j] == Types) { INTEGER(Counts)[8] +=1;
                            } else {INTEGER(Counts)[6] +=1; }
                        } else {
                            if (X[j] == Types) { INTEGER(Counts)[7] +=1;
                            } else {
                                if (X[j] == F[j]) { INTEGER(Counts)[4] +=1;
                                } else {INTEGER(Counts)[5] +=1; }
                                
                            }
                        }
                    }
                }
            }
        }
        
        UNPROTECT(1);
        return Counts;
    }
}
