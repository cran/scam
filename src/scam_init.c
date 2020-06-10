/* Symbol registration initialization: original provided by Brian Ripley.
   Anything called from R should be registered here (and declared in scam.h).
   (See also NAMESPACE:1)
 */ 
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "scam.h"


R_CMethodDef CEntries[] = { 
    {"rwMatrix",(DL_FUNC)&rwMatrix,8},
    {NULL, NULL, 0}
};

void R_init_scam(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);   
}
