#include <stdio.h>
#include <stdlib.h>
#include "perf.h"
#include "cblas.h"


int main(){
  double performance;
  perf_t start,stop;

  double a,b,c;
  a = 1;
  b = 1;
  long flop = 1;

  // Executions a vide, flush potentiel, ...

  // Performance d'une addition scalaire
  perf(&start);
  c = a + b;
  perf(&stop);

  // Verification
  printf("%lf = %lf + %lf\n", c, a, b);

  // Performance
  perf_diff(&start, &stop);
  performance = perf_mflops(&stop, flop);
  printf("Mflop/s : %lf \n", performance);

  return 0;
}
