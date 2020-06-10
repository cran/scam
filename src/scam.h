/* main method routines */

/* See http://developer.r-project.org/blosxom.cgi/R-devel/2019/08/29#n2019-08-29
   for what USE_FC_LEN_T is doing and for why see
   https://developer.r-project.org/Blog/public/2019/05/15/gfortran-issues-with-lapack/index.html 
   
   In a nutshell, the mechanism used to call BLAS/LAPACK from C (by everyone, not just R) is not 
   technically supported by the Fortran standard. Fortran needs to know how long strings are (they 
   are not null terminated) so the lengths are passed as hidden extra function arguments. BLAS/LAPACK
   only ever uses single character strings, so it never needs to access the string lengths and it is 
   then no problem that they are missing (they are at the end of the call arguments), so they are simply 
   not passed in the C call. This was no problme until Gfortran decided to optimize the process of calling 
   a function with the same argument list as the calling function. Basically it passed the call stack of 
   the calling function to the called function assuming that it contained the string lengths - as it 
   didn't this caused stack corruption. 

   The solution is to pass the unit string lengths explicitly using FCONE defined in Rconfig.h if 
   USE_FC_LEN_T is defined. This mechanism is needed since it is compiler specific what type is used 
   to pass the string lengths (what happens then if BLAS/LAPACK and R are compiled using different 
   extra argument types is unclear to me, but no problems of this sort are currently known in any case 
   to get an actual problem the LAPACK/BLAS compiler would have to be using a different number of bytes 
   to the R compiler). 

   In practice when calling BLAS/LAPACK macro FCONE has to be added to the end of the call as
   many times as there are character arguments to the call. mat.c has many examples.

*/

#define USE_FC_LEN_T
#include <Rinternals.h>
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
/* If we are compiling with a version of R before FCONE and the explicit supplying of extra arguments 
   was introduced, then FCONE has to be defined */ 
#ifndef FCONE
#define FCONE
#endif

/* Most compilers with openMP support supply 
   a pre-defined compiler macro _OPENMP. Following 
   facilitates selective turning off (by testing value 
   or defining multiple versions OPENMP_ON1, OPENMP_ON2...)  */

#if defined _OPENMP
#define OPENMP_ON 1 
#endif
/* ... note also that there is no actual *need* to protect #pragmas with 
  #ifdef OPENMP_ON, since C ignores undefined pragmas, but failing 
  to do so may produce alot of compilation warnings if openMP is not supported. 
  In contrast functions from omp.h must be protected, and there is 
  non-avoidable use of these in the mgcv code. */

//#define OMP_REPORT // define to have all routines using omp report on start and end.

/* sed -i 's/old-text/new-text/g' *.c
   is quite useful!!
*/

/* For safe memory handling from R... */
#define CALLOC R_chk_calloc
#define FREE R_chk_free
/* BUT, this can mess up valgrinding for memory error checking - problems are 
   sometimes missed because standard allocation is being circumvented. Then errors can 
   corrupt R memory management without detection and trigger nothing until R
   messes up internally because of corruption, which then makes it look as if
   R is generating the problem. Hence better to reset for checking. Also sizing
   errors in .C often generate no obvious valgrind error.*/
//#define CALLOC calloc
//#define FREE free

/* void *R_chk_calloc1(size_t nmemb,size_t size);

/* service routines */

void rwMatrix(int *stop,int *row,double *w,double *X,int *n,int *p,int *trans,double *work);

