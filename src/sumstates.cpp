/*
 * File sumstates.cpp
 *
 * Part of the `epiPOMS' R package
 *
 * Calculation of total number of times than an individual belongs to a specific state over the
 * entire time period
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    
    SEXP sumstates(SEXP B1, SEXP I1, SEXP Types1){
        int *B;
        int  j, i, I, Types;
        SEXP infected_states;
        
        B1 = coerceVector(B1, INTSXP);
        Types1 = coerceVector(Types1, INTSXP);
        I1 = coerceVector(I1, INTSXP);
        
        B = INTEGER(B1);
        Types = INTEGER(Types1)[0];
        I = INTEGER(I1)[0];
        
        PROTECT(infected_states = allocVector(INTSXP, Types));
        
        //  set to zero
        for (j=0 ; j<Types ; j++) {
            INTEGER(infected_states)[j] = 0;
        }
        
        //  for all individuals find states
        for (i =0 ; i<I ; i++) {
            for (j=0 ; j<Types ; j++) {
                if ( B[i]==j) {
                    INTEGER(infected_states)[j] +=1;
                }
            }
        }
        
        
        UNPROTECT(1);
        return infected_states;
    }
}
