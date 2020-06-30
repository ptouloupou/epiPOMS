#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP condlike(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP gradient(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP logpost(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP statetrans(SEXP, SEXP, SEXP);
extern SEXP sumstates(SEXP, SEXP, SEXP);
extern SEXP transprob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP transprobnotNA(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP updateobspars(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"condlike",       (DL_FUNC) &condlike,        7},
    {"gradient",       (DL_FUNC) &gradient,       14},
    {"logpost",        (DL_FUNC) &logpost,        12},
    {"statetrans",     (DL_FUNC) &statetrans,      3},
    {"sumstates",      (DL_FUNC) &sumstates,       3},
    {"transprob",      (DL_FUNC) &transprob,       9},
    {"transprobnotNA", (DL_FUNC) &transprobnotNA, 10},
    {"updateobspars",  (DL_FUNC) &updateobspars,   5},
    {NULL, NULL, 0}
};

void R_init_epiPOMS(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
