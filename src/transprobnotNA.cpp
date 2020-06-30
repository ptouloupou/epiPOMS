/*
 * File transprobnotNA.cpp
 *
 * Part of the `epiPOMS' R package
 *
 * Calculation of the transition probabilities in the case of a missing individual required by
 * iFFBS
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    
    SEXP transprobnotNA(SEXP a1, SEXP b1, SEXP m1, SEXP delta1, SEXP gamma1, SEXP I1, SEXP B1, SEXP Types1, SEXP zerosint, SEXP counts1){
        int *B, *counts, *infected_states;
        double  delta, gamma;
        double sumInfrate,  sumdeltaInfrateneqj1,  probhelp1;
        int    j, j1, j2, j3, i,  I, Types;
        double *a, *b, *m;
        SEXP flag, logpost, Trans, Rec, Infec, Infectplusone;
        
        B1 = coerceVector(B1, INTSXP);
        I1 = coerceVector(I1, INTSXP);
        counts1 = coerceVector(counts1, INTSXP);
        Types1 = coerceVector(Types1, INTSXP);
        a1 = coerceVector(a1, REALSXP);
        b1 = coerceVector(b1, REALSXP);
        m1 = coerceVector(m1, REALSXP);
        delta1 = coerceVector(delta1, REALSXP);
        gamma1 = coerceVector(gamma1, REALSXP);
        zerosint = coerceVector(zerosint, INTSXP);
        
        
        B = INTEGER(B1);
        counts = INTEGER(counts1);
        I = INTEGER(I1)[0];
        Types = INTEGER(Types1)[0];
        delta = REAL(delta1)[0];
        gamma = REAL(gamma1)[0];
        infected_states = INTEGER(zerosint);
        a = REAL(a1);
        b = REAL(b1);
        m = REAL(m1);
        
        PROTECT(logpost = allocVector(REALSXP, Types + Types ));
        PROTECT(Trans = allocVector(REALSXP, (pow(Types,2)) ));
        PROTECT(Rec = allocVector(REALSXP, Types));
        PROTECT(Infec = allocVector(REALSXP, Types));
        PROTECT(flag = allocVector(REALSXP, Types));
        PROTECT(Infectplusone = allocVector(INTSXP, Types));
        
        //  set to zero or one
        for (j=0 ; j<Types ; j++) {
            infected_states[j] = 0;
            REAL(logpost)[j] = 1.0;
            REAL(flag)[j] = 1.0;
            INTEGER(Infectplusone)[j] = 0;
        }
        
        //  for all individuals find states
        for (i =0 ; i<I ; i++) {
            for (j=0 ; j<Types ; j++) {
                if ( B[i]==j) {
                    infected_states[j] +=1;
                }
            }
        }
        
        for (j=0 ; j< (pow(Types,2)) ; j++) {
            counts[j] = 0;
        }
        
        //  calculate the counts
        for (i=0 ; i<I ; i++) {
            //  double for loop for types
            for (j1=0 ; j1<Types ; j1++) {
                for (j2=0 ; j2<Types ; j2++) {
                    if ((B[i]==j1) && (B[I+i]==j2)) {
                        counts[j1*Types+j2] += 1;
                    }
                }
            }
        }
        
        for (j3=0 ; j3<Types ; j3++) {
            
            for (j=0 ; j<Types ; j++) {
                if (j3!=j) {
                    INTEGER(Infectplusone)[j] = infected_states[j] ;
                } else {
                    INTEGER(Infectplusone)[j] = infected_states[j] + 1 ;
                }
            }
            
            
            for (j=0 ; j< (pow(Types,2)) ; j++) {
                REAL(Trans)[j] = 0.0;
            }
            
            sumInfrate = 0.0;
            
            for (j=1 ; j<Types ; j++) {
                
                REAL(Infec)[j] =  a[j] + INTEGER(Infectplusone)[j]*b[j]*gamma;
                REAL(Rec)[j] = m[j];
                
                sumInfrate += REAL(Infec)[j];
                
                //Rprintf("\n values is %.5f", DlogInfecb[j]);
                // Rprintf("\n values is %d", infected_states[j]);
            }
            
            //  Calculate Transision matrix
            REAL(Trans)[0] = exp(-sumInfrate);
            for (j2=1 ; j2<Types ; j2++) {
                probhelp1 = (1-exp( - sumInfrate)) / ( sumInfrate);
                REAL(Trans)[j2] = (REAL(Infec)[j2] )  * probhelp1;
            }
            
            for (j1=1 ; j1<Types ; j1++) {
                sumdeltaInfrateneqj1 = 0.0;
                
                for (j2=1 ; j2<Types ; j2++) {
                    if (j2!=j1) {
                        sumdeltaInfrateneqj1 += delta * REAL(Infec)[j2];
                    }
                }
                
                REAL(Trans)[j1*Types+j1] = exp( -m[j1] - sumdeltaInfrateneqj1);
                probhelp1 = (1-exp( -m[j1] - sumdeltaInfrateneqj1)) / ( m[j1] + sumdeltaInfrateneqj1);
                
                REAL(Trans)[j1*Types] = (m[j1] ) * probhelp1;
                
                for (j2=1 ; j2<Types ; j2++) {
                    if (j2!=j1) {
                        REAL(Trans)[j1*Types+j2] = ( delta * REAL(Infec)[j2])  * probhelp1;       
                    }
                }
                
            }
            
            
            for (j=0 ; j< (pow(Types,2)) ; j++) {
                if (REAL(Trans)[j] < 0.0){ REAL(flag)[j3] = 0.0;}
            }
            
            for (j=0 ; j< (pow(Types,2)) ; j++) {
                REAL(logpost)[j3] = REAL(logpost)[j3] * pow( REAL(Trans)[j] , counts[j]); 
            }
        }
        
        for (j=Types ; j< (2*Types) ; j++) {
            REAL(logpost)[j] = REAL(flag)[j-Types] ; 
        }
        
        UNPROTECT(6);
        return logpost;
    }        
}      
