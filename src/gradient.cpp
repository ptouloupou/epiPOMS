/*
 * File gradient.cpp
 *
 * Part of the `epiPOMS' R package
 *
 * Calculation of the gradient expressions for parameters updated with HMC
 */


#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

extern "C" {
    
    SEXP gradient(SEXP a1, SEXP b1, SEXP m1, SEXP delta1,  SEXP gamma1, SEXP I1, SEXP T1, SEXP P1, SEXP B1, SEXP Types1, SEXP zerosint, SEXP counts1, SEXP alphas1, SEXP Last1){
        int *B, *counts, *infected_states;
        double sumInfrate, sumInfrateneqj1, sumdeltaInfrateneqj1, sumdeltaInfrateneqj2, probhelp1, probhelp2;
        int sumcountsneqj1, sumcountsneqj2, sumcounts0;
        double delta, gamma;
        int  tmp1, j, jc, j1, j2, j3, i, t, I, T, P, Types, Last;
        double *a, *b, *m, *alphas;
        SEXP theta, Infrate, DInfratea, DInfrateb, DInfrategamma, Recrate, DRecrate, Ddelta;
        
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
        alphas1 = coerceVector(alphas1, REALSXP);
        Last1 = coerceVector(Last1, INTSXP);
        
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
        alphas = REAL(alphas1);
        Last = INTEGER(Last1)[0];
        
        PROTECT(theta = allocVector(REALSXP, 3*(Types-1) + 2));
        PROTECT(Infrate = allocVector(REALSXP, Types));
        PROTECT(DInfratea = allocVector(REALSXP, Types));
        PROTECT(DInfrateb = allocVector(REALSXP, Types));
        PROTECT(DInfrategamma = allocVector(REALSXP, Types));
        PROTECT(Recrate = allocVector(REALSXP, Types));
        PROTECT(DRecrate = allocVector(REALSXP, Types));
        PROTECT(Ddelta = allocVector(REALSXP, 1));
        
        
        for (j=0 ; j<(3*(Types-1)+2) ; j++) {
            REAL(theta)[j] = 0.0;
        }
        
        
        if (P == 0)
        {
            for (t=1 ; t<T ; t++) {
                //  set to zero
                for (j=0 ; j<Types ; j++) {
                    infected_states[j] = 0;
                }
                //  for all individuals find states
                for (i =0 ; i<I ; i++) {
                    for (j=0 ; j<Types ; j++) {
                        if ( B[(t-1)*I+i]==j) {
                            infected_states[j] +=1;
                        }
                    }
                }
                sumInfrate = 0.0;
                
                for (j=1 ; j<Types ; j++) {
                    
                    REAL(Infrate)[j] =  a[j] + infected_states[j]*b[j];
                    REAL(DInfratea)[j] = a[j];
                    REAL(DInfrateb)[j] = infected_states[j]*b[j];
                    REAL(Recrate)[j] = m[j];
                    REAL(DRecrate)[j] = m[j];
                    REAL(Ddelta)[0] = delta;
                    
                    sumInfrate += REAL(Infrate)[j];
                    
                }
                
                for (j=0 ; j< (pow(Types,2)) ; j++) {
                    counts[j] = 0;
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
                
                //  Gradient part A
                
                for (j1=1 ; j1<Types ; j1++) {
                    //  A1
                    sumInfrateneqj1 = 0.0;
                    sumdeltaInfrateneqj1 = 0.0;
                    sumcountsneqj1 = 0;
                    tmp1=0;
                    
                    for (j2=1 ; j2<Types ; j2++) {
                        sumdeltaInfrateneqj2 = 0.0;
                        sumcountsneqj2 = 0;
                        
                        if (j2!=j1) {
                            sumInfrateneqj1 += REAL(Infrate)[j2];
                            sumdeltaInfrateneqj1 += delta * REAL(Infrate)[j2];
                            sumcountsneqj1 += counts[j1*Types+j2];
                            tmp1 += counts[j2*Types+j1];
                            
                            for (j3=1 ; j3<Types ; j3++) {
                                if (j3!=j2) {
                                    sumdeltaInfrateneqj2 += delta * REAL(Infrate)[j3];
                                    sumcountsneqj2 += counts[j2*Types+j3];
                                }
                            }
                            
                            // derivate of alpha's, beta's
                            REAL(theta)[j1-1] += 0 - counts[j2*Types+j2] * delta * REAL(DInfratea)[j1];
                            REAL(theta)[(Types-1)+j1-1] += 0 - counts[j2*Types+j2] * delta * REAL(DInfrateb)[j1];
                            
                            probhelp1 = exp(-m[j2] -  sumdeltaInfrateneqj2) / (1 - exp(-m[j2] -  sumdeltaInfrateneqj2));
                            probhelp2 = 1 / ( m[j2] +  sumdeltaInfrateneqj2 );
                            REAL(theta)[j1-1] += (sumcountsneqj2 + counts[j2*Types]) * delta * (REAL(DInfratea)[j1]) * ( probhelp1 - probhelp2);
                            REAL(theta)[(Types-1)+j1-1] += (sumcountsneqj2 + counts[j2*Types]) * delta * (REAL(DInfrateb)[j1]) * ( probhelp1 - probhelp2);
                            
                            // derivate of delta
                            probhelp1 = 1/ (delta * REAL(Infrate)[j2]);
                            REAL(theta)[3*(Types-1)] += counts[j1*Types+j2] *  REAL(Infrate)[j2] * REAL(Ddelta)[0] * probhelp1;
                        }
                    }
                    
                    // derivate of mu's
                    probhelp1 = 1 / (REAL(Recrate)[j1]);
                    REAL(theta)[2*(Types-1)+j1-1] += counts[j1*Types] * REAL(DRecrate)[j1] * probhelp1 - counts[j1*Types + j1] * REAL(DRecrate)[j1];
                    probhelp1 = exp(-m[j1] -  sumdeltaInfrateneqj1) / (1 - exp(-m[j1] -  sumdeltaInfrateneqj1));
                    probhelp2 = 1 / ( m[j1] +  sumdeltaInfrateneqj1);
                    REAL(theta)[2*(Types-1)+j1-1] += (sumcountsneqj1 + counts[j1*Types]) * (REAL(DRecrate)[j1]) * ( probhelp1 - probhelp2);
                    
                    // derivate of delta
                    REAL(theta)[3*(Types-1)] += 0.0 - counts[j1*Types + j1] * REAL(Ddelta)[0] * sumInfrateneqj1;
                    probhelp1 = exp(-m[j1] -  sumdeltaInfrateneqj1) / (1 - exp(-m[j1] -  sumdeltaInfrateneqj1));
                    probhelp2 = 1 / ( m[j1] +  sumdeltaInfrateneqj1);
                    REAL(theta)[3*(Types-1)] += (sumcountsneqj1 + counts[j1*Types]) * ( REAL(Ddelta)[0] ) *(sumInfrateneqj1) * (  probhelp1 - probhelp2);
                    
                    // derivate of alpha's, beta's
                    sumcounts0 = 0;
                    for (jc=1 ; jc<Types ; jc++) {
                        sumcounts0 += counts[jc];
                    }
                    
                    probhelp1 = 1/ (delta * REAL(Infrate)[j1]);
                    REAL(theta)[j1-1] += tmp1 * delta * REAL(DInfratea)[j1]  * probhelp1;
                    REAL(theta)[(Types-1)+j1-1] += tmp1 * delta * REAL(DInfrateb)[j1]  * probhelp1;
                    
                    probhelp1 = exp(-sumInfrate) / (1 - exp(-sumInfrate));
                    probhelp2 = 1 / sumInfrate;
                    REAL(theta)[j1-1] += sumcounts0 * REAL(DInfratea)[j1] * ( probhelp1 - probhelp2) - counts[0] * REAL(DInfratea)[j1] ;
                    REAL(theta)[(Types-1)+j1-1] += sumcounts0 * REAL(DInfrateb)[j1] * ( probhelp1 - probhelp2) - counts[0] * REAL(DInfrateb)[j1] ;
                    
                    probhelp1 = 1 / (REAL(Infrate)[j1]);
                    REAL(theta)[j1-1] +=  counts[j1]* REAL(DInfratea)[j1] * (probhelp1);
                    REAL(theta)[(Types-1)+j1-1] += counts[j1]* REAL(DInfrateb)[j1] * (probhelp1);
                    
                    
                }
                
                
                
            } // close time points
            
            
        } else {
            for (t=1 ; t<T ; t++) {
                //  set to zero
                for (j=0 ; j<Types ; j++) {
                    infected_states[j] = 0;
                }
                //  for all individuals find states
                for (i =0 ; i<I ; i++) {
                    for (j=0 ; j<Types ; j++) {
                        if ( B[(t-1)*I+i]==j) {
                            infected_states[j] +=1;
                        }
                    }
                }
                sumInfrate = 0.0;
                
                for (j=1 ; j<Types ; j++) {
                    
                    REAL(Infrate)[j] =  a[j] + infected_states[j]*b[j]*gamma;
                    REAL(DInfratea)[j] = a[j];
                    REAL(DInfrateb)[j] = gamma*infected_states[j]*b[j];
                    REAL(DInfrategamma)[j] = gamma*infected_states[j]*b[j];
                    REAL(Recrate)[j] = m[j];
                    REAL(DRecrate)[j] = m[j];
                    REAL(Ddelta)[0] = delta;
                    
                    
                    sumInfrate += REAL(Infrate)[j];
                    
                }
                
                for (j=0 ; j< (pow(Types,2)) ; j++) {
                    counts[j] = 0;
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
                
                //  Gradient part A
                
                for (j1=1 ; j1<Types ; j1++) {
                    //  A1
                    sumInfrateneqj1 = 0.0;
                    sumdeltaInfrateneqj1 = 0.0;
                    sumcountsneqj1 = 0;
                    tmp1=0;
                    
                    for (j2=1 ; j2<Types ; j2++) {
                        sumdeltaInfrateneqj2 = 0.0;
                        sumcountsneqj2 = 0;
                        
                        if (j2!=j1) {
                            sumInfrateneqj1 += REAL(Infrate)[j2];
                            sumdeltaInfrateneqj1 += delta * REAL(Infrate)[j2];
                            sumcountsneqj1 += counts[j1*Types+j2];
                            tmp1 += counts[j2*Types+j1];
                            
                            for (j3=1 ; j3<Types ; j3++) {
                                if (j3!=j2) {
                                    sumdeltaInfrateneqj2 += delta * REAL(Infrate)[j3];
                                    sumcountsneqj2 += counts[j2*Types+j3];
                                }
                            }
                            
                            // derivate of alpha's, beta's and gamma
                            REAL(theta)[j1-1] += 0 - counts[j2*Types+j2] * delta * REAL(DInfratea)[j1];
                            REAL(theta)[(Types-1)+j1-1] += 0 - counts[j2*Types+j2] * delta * REAL(DInfrateb)[j1];
                            REAL(theta)[3*(Types-1) + 1] += 0 - counts[j2*Types+j2] * delta * REAL(DInfrategamma)[j1];
                            
                            probhelp1 = exp(-m[j2] -  sumdeltaInfrateneqj2) / (1 - exp(-m[j2] -  sumdeltaInfrateneqj2));
                            probhelp2 = 1 / ( m[j2] +  sumdeltaInfrateneqj2 );
                            REAL(theta)[j1-1] += (sumcountsneqj2 + counts[j2*Types]) * delta * (REAL(DInfratea)[j1]) * ( probhelp1 - probhelp2);
                            REAL(theta)[(Types-1)+j1-1] += (sumcountsneqj2 + counts[j2*Types]) * delta * (REAL(DInfrateb)[j1]) * ( probhelp1 - probhelp2);
                            REAL(theta)[3*(Types-1) + 1] += (sumcountsneqj2 + counts[j2*Types]) * delta * (REAL(DInfrategamma)[j1]) * ( probhelp1 - probhelp2);
                            
                            // derivate of delta
                            probhelp1 = 1/ (delta * REAL(Infrate)[j2]);
                            REAL(theta)[3*(Types-1)] += counts[j1*Types+j2] *  REAL(Infrate)[j2] * REAL(Ddelta)[0] * probhelp1;
                        }
                    }
                    
                    // derivate of mu's
                    probhelp1 = 1 / (REAL(Recrate)[j1]);
                    REAL(theta)[2*(Types-1)+j1-1] += counts[j1*Types] * REAL(DRecrate)[j1] * probhelp1 - counts[j1*Types + j1] * REAL(DRecrate)[j1];
                    probhelp1 = exp(-m[j1] -  sumdeltaInfrateneqj1) / (1 - exp(-m[j1] -  sumdeltaInfrateneqj1));
                    probhelp2 = 1 / ( m[j1] +  sumdeltaInfrateneqj1);
                    REAL(theta)[2*(Types-1)+j1-1] += (sumcountsneqj1 + counts[j1*Types]) * (REAL(DRecrate)[j1]) * ( probhelp1 - probhelp2);
                    
                    // derivate of delta
                    REAL(theta)[3*(Types-1)] += 0.0 - counts[j1*Types + j1] * REAL(Ddelta)[0] * sumInfrateneqj1;
                    probhelp1 = exp(-m[j1] -  sumdeltaInfrateneqj1) / (1 - exp(-m[j1] -  sumdeltaInfrateneqj1));
                    probhelp2 = 1 / ( m[j1] +  sumdeltaInfrateneqj1);
                    REAL(theta)[3*(Types-1)] += (sumcountsneqj1 + counts[j1*Types]) * ( REAL(Ddelta)[0] ) *(sumInfrateneqj1) * (  probhelp1 - probhelp2);
                    
                    // derivate of alpha's, beta's and gamma
                    sumcounts0 = 0;
                    for (jc=1 ; jc<Types ; jc++) {
                        sumcounts0 += counts[jc];
                    }
                    
                    probhelp1 = 1/ (delta * REAL(Infrate)[j1]);
                    REAL(theta)[j1-1] += tmp1 * delta * REAL(DInfratea)[j1]  * probhelp1;
                    REAL(theta)[(Types-1)+j1-1] += tmp1 * delta * REAL(DInfrateb)[j1]  * probhelp1;
                    REAL(theta)[3*(Types-1) + 1] += tmp1 * delta * REAL(DInfrategamma)[j1]  * probhelp1;
                    
                    probhelp1 = exp(-sumInfrate) / (1 - exp(-sumInfrate));
                    probhelp2 = 1 / sumInfrate;
                    REAL(theta)[j1-1] += sumcounts0 * REAL(DInfratea)[j1] * ( probhelp1 - probhelp2) - counts[0] * REAL(DInfratea)[j1] ;
                    REAL(theta)[(Types-1)+j1-1] += sumcounts0 * REAL(DInfrateb)[j1] * ( probhelp1 - probhelp2) - counts[0] * REAL(DInfrateb)[j1] ;
                    REAL(theta)[3*(Types-1) + 1] += sumcounts0 * REAL(DInfrategamma)[j1] * ( probhelp1 - probhelp2) - counts[0] * REAL(DInfrategamma)[j1] ;
                    
                    
                    probhelp1 = 1 / (REAL(Infrate)[j1]);
                    REAL(theta)[j1-1] +=  counts[j1]* REAL(DInfratea)[j1] * (probhelp1);
                    REAL(theta)[(Types-1)+j1-1] += counts[j1]* REAL(DInfrateb)[j1] * (probhelp1);
                    REAL(theta)[3*(Types-1) + 1] += counts[j1]* REAL(DInfrategamma)[j1] * (probhelp1);
                    
                }
            } //close time points
        } //close else
        
        if (Last == 1) {
            REAL(theta)[3*(Types-1)] += 1 - alphas[3]*delta;
            REAL(theta)[3*(Types-1) + 1] += 1 - alphas[4]*gamma;
            
            for (j=1 ; j<Types ; j++) {
                REAL(theta)[j-1] +=  1 - alphas[0] * a[j];
                REAL(theta)[(Types-1)+j-1] +=  1 - alphas[1] * b[j];
                REAL(theta)[2*(Types-1)+j-1] += 1 - alphas[2]  * m[j];
            }
        }
        
        UNPROTECT(8);
        return theta;
    }
}
