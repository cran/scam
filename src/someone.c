/* Copyright (C) 2008-2014 Simon N. Wood  simon.wood@r-project.org

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
(www.gnu.org/copyleft/gpl.html)

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
USA. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rconfig.h>
#include "scam.h"



/*******************************************************/
/** Fast re-weighting routines                         */
/*******************************************************/

void rwMatrix(int *stop,int *row,double *w,double *X,int *n,int *p,int *trans,double *work) {
/* Function to recombine rows of n by p matrix X (column ordered).
   ith row of X' is made up of row[stop[i-1]+1...stop[i]], weighted by 
   w[stop[i-1]+1...stop[i]]. stop[-1]=-1 by convention.
   stop is an n vector.     
   
   If (trans==0) the operation on a column x is x'[i] += w[row[j]] * X[row[j]] over the 
   j from stop[i-1]+1 to stop[i]. Otherwise the tranposed operation 
   x'[row[j]] += w[row[j]] * x[i] is used with the same j range. x' zero at outset.

   work is same dimension as X

   See rwMatrix in bam.r for call from R. 
*/
  ptrdiff_t i,j,jump,start=0,end,off;
  double *X1p,*Xp,weight,*Xpe,*X1;
  /* create storage for output matrix, cleared to zero */
  X1 = work; 
  jump = *n;
  off = *n * (ptrdiff_t) *p; 
  for (X1p=X1,Xpe=X1p+off;X1p<Xpe;X1p++) *X1p = 0.0;
  for (i=0;i<*n;i++) { /* loop through rows of output X1 */
    end = stop[i]+1;
    for (j=start;j<end;j++) { /* loop through the input rows */
      if (*trans) {
        X1p = X1 + row[j];
        Xp = X + i;
      } else { 
        X1p = X1 + i;    /* pointer to start of row i of output */
        Xp = X + row[j]; /* pointer to start of source row */
      }
      weight = w[j];   
      for (Xpe=Xp+off;Xp<Xpe;Xp+=jump,X1p+=jump) *X1p += weight * *Xp;
    }
    start = end;
  }
  /* copy output to input for return...*/
  for (Xp=X,X1p=X1,Xpe=Xp+off;Xp<Xpe;Xp++,X1p++) *Xp = *X1p;
}

/* Example code for rwMatrix in R....
   n <- 10;p<-5
   X <- matrix(runif(n*p),n,p)
   ## create transform to take AR1(rho) to independence...
   stop <- c(1:(n-1)*2,2*n-1)
   row <- rep(1:n,rep(2,n))[-1]
   rho <- .7;ld <- 1/sqrt(1-rho^2);sd <- -rho*ld
   w <- c(rep(c(ld,sd),n-1),1)
   mgcv:::rwMatrix(stop,row,w,X)
   
*/

