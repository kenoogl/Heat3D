
#include <stdio.h>
#include <stdlib.h>
#include "stretch.h"

int main(int argc, char * argv[]) {

  // Check command line arguments
  if ( argc != 6) {
    printf("Usage:\n");
    printf("\t$ ./test size x1 x2 sp1 sp2\n");
    return -1;
  }

  int size;
  double sp1, sp2;

  double x1;
  double x2;

  size = atoi(argv[1]);
  x1  = atof(argv[2]);
  x2  = atof(argv[3]);
  sp1 = atof(argv[4]);
  sp2 = atof(argv[5]);

  printf("start = %lf\n", x1); // 開始点(x1)
  printf("end   = %lf\n", x2); // 終了点(x2)
  printf("size = %d\n", size); // 分割ノード点
  printf("sp1  = %lf\n", sp1); // 開始点(x1)側の格子点間隔
  printf("sp2  = %lf\n", sp2); // 終了点(x2)側の格子点間隔

  REAL_TYPE* x = new REAL_TYPE[size];

  Stretch s(size, x1, x2, sp1, sp2);

  if ( !s.distribution(x) )
  {
    printf("Not converged.\n");
  };

  delete [] x;

  return 0;
}
