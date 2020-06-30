/*
 * File logpost.cpp
 *
 * Part of the `epiPOMS' R package
 *
 * Calculation of the log posteriors required by HMC
 */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    
    SEXP logpost(SEXP a1, SEXP b1, SEXP m1, SEXP delta1, SEXP gamma1, SEXP I1, SEXP T1, SEXP P1, SEXP B1, SEXP Types1, SEXP zerosint, SEXP counts1){
        int *B, *counts, *infected_states;
        double  delta, gamma;
        double sumInfrate, sumdeltaInfrateneqj1, sumInfrateneqj1, probhelp1;
        int  j, j1, j2, i, t, I, T, P, Types;
        double *a, *b, *m;
        SEXP logpost, Trans, Rec, Infec;
        
        B1 = coerceVector(B1, INTSXP);
        I1 = coerceVector(I1, INTSXP);
        T1 = coerceVector(T1, INTSXP);
        P1 = coerceVector(P1, INTSXP);
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
        T = INTEGER(T1)[0];
        P = INTEGER(P1)[0];
        Types = INTEGER(Types1)[0];
        delta = REAL(delta1)[0];
        gamma = REAL(gamma1)[0];
        infected_states = INTEGER(zerosint);
        a = REAL(a1);
        b = REAL(b1);
        m = REAL(m1);
        
        PROTECT(logpost = allocVector(REALSXP, 1));
        PROTECT(Trans = allocVector(REALSXP, (pow(Types,2)) ));
        PROTECT(Rec = allocVector(REALSXP, Types));
        PROTECT(Infec = allocVector(REALSXP, Types));
        
        REAL(logpost)[0] = 0.0;
        
        if (P == 0)
        {
            for (t=1 ; t<T ; t++) {
                //  set to zero
                for (j=0 ; j<Types ; j++) {
                    infected_states[j] = 0;
                }
                
                for (j=0 ; j< (pow(Types,2)) ; j++) {
                    counts[j] = 0;
                    REAL(Trans)[j] = 0.0;
                }
                
                //  for all individuals find states
                for (i =0 ; i<I ; i++) {
                    for (j=0 ; j<Types ; j++) {
                        if ( B[(t-1)*I+i]==j) {
                            infected_states[j] +=1;
                        }
                    }
                }
                
                //  calculate the counts
                for (i=0 ; i<I ; i++) {
                    //  double for loop for types
                    for (j1=0 ; j1<Types ; j1++) {
                        for (j2=0 ; j2<Types ; j2++) {
                            if ((B[(t-1)*I+i]==j1) && (B[t*I+i]==j2)) {
                                counts[j1*Types+j2] += 1;
                            }
                        }
                    }
                }
                
                sumInfrate = 0.0;
                
                for (j=1 ; j<Types ; j++) {
                    
                    REAL(Infec)[j] =  a[j] + infected_states[j]*b[j];
                    REAL(Rec)[j] = m[j];
                    
                    sumInfrate += REAL(Infec)[j];
                    
                }
                
                //  Calculate Transision matrix
                REAL(Trans)[0] = exp(-sumInfrate);
                for (j2=1 ; j2<Types ; j2++) {
                    probhelp1 = (1-exp( - sumInfrate)) / ( sumInfrate );
                    REAL(Trans)[j2] = (REAL(Infec)[j2] ) * probhelp1;
                }
                
                for (j1=1 ; j1<Types ; j1++) {
                    sumdeltaInfrateneqj1 = 0.0;
                    
                    for (j2=1 ; j2<Types ; j2++) {
                        if (j2!=j1) {
                            sumdeltaInfrateneqj1 += delta * REAL(Infec)[j2];
                        }
                    }
                    
                    REAL(Trans)[j1*Types+j1] = exp( -m[j1] - sumdeltaInfrateneqj1);
                    probhelp1 = (1-exp( -m[j1] - sumdeltaInfrateneqj1)) / ( sumdeltaInfrateneqj1 +m[j1] );
                    
                    REAL(Trans)[j1*Types] = (m[j1] ) * probhelp1;
                    
                    
                    for (j2=1 ; j2<Types ; j2++) {
                        if (j2!=j1) {
                            REAL(Trans)[j1*Types+j2] = ( delta * REAL(Infec)[j2] ) * probhelp1;
                            
                        }
                    }
                    
                }
                
                
                for (j=0 ; j< (pow(Types,2)) ; j++) {
                    if (counts[j]>0) {
                        REAL(logpost)[0] += counts[j]*log(REAL(Trans)[j]);
                    }
                }
            } // close time points
            
        } else {
            for (t=1 ; t<T ; t++) {
                //  set to zero
                for (j=0 ; j<Types ; j++) {
                    infected_states[j] = 0;
                }
                
                for (j=0 ; j< (pow(Types,2)) ; j++) {
                    counts[j] = 0;
                    REAL(Trans)[j] = 0.0;
                }
                
                //  for all individuals find states
                for (i =0 ; i<I ; i++) {
                    for (j=0 ; j<Types ; j++) {
                        if ( B[(t-1)*I+i]==j) {
                            infected_states[j] +=1;
                        }
                    }
                }
                
                //  calculate the counts
                for (i=0 ; i<I ; i++) {
                    //  double for loop for types
                    for (j1=0 ; j1<Types ; j1++) {
                        for (j2=0 ; j2<Types ; j2++) {
                            if ((B[(t-1)*I+i]==j1) && (B[t*I+i]==j2)) {
                                counts[j1*Types+j2] += 1;
                            }
                        }
                    }
                }
                
                sumInfrate = 0.0; 
                
                for (j=1 ; j<Types ; j++) {
                    
                    REAL(Infec)[j] =  a[j] + gamma*infected_states[j]*b[j];
                    REAL(Rec)[j] = m[j];
                    
                    sumInfrate += REAL(Infec)[j];
                    
                }
                
                //  Calculate Transision matrix
                REAL(Trans)[0] = exp(-sumInfrate); 
                for (j2=1 ; j2<Types ; j2++) {
                    probhelp1 = (1-exp( - sumInfrate)) / ( sumInfrate);
                    REAL(Trans)[j2] = (REAL(Infec)[j2] ) * probhelp1; 
                }
                
                for (j1=1 ; j1<Types ; j1++) {
                    sumdeltaInfrateneqj1 = 0.0;
                    sumInfrateneqj1 = 0.0;
                    
                    for (j2=1 ; j2<Types ; j2++) {
                        if (j2!=j1) {
                            sumInfrateneqj1 += REAL(Infec)[j2];
                            sumdeltaInfrateneqj1 += delta * REAL(Infec)[j2];
                        }
                    }
                    
                    REAL(Trans)[j1*Types+j1] = exp( -m[j1] - sumdeltaInfrateneqj1);
                    probhelp1 = (1-exp( -m[j1] - sumdeltaInfrateneqj1)) / ( m[j1] + sumdeltaInfrateneqj1 );
                    
                    REAL(Trans)[j1*Types] = (m[j1] ) * probhelp1;
                    
                    
                    for (j2=1 ; j2<Types ; j2++) {
                        if (j2!=j1) {
                            REAL(Trans)[j1*Types+j2] = ( delta * REAL(Infec)[j2] ) * probhelp1;       
                        }
                    }
                    
                }
                
                for (j=0 ; j< (pow(Types,2)) ; j++) {
                    
                    if (counts[j]>0) {
                        REAL(logpost)[0] += counts[j]*log(REAL(Trans)[j]);
                    }
                    
                }
            } // close time points
        } //close if 
        
        
        UNPROTECT(4);
        return logpost;
    }        
}      
