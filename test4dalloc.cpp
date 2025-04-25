#include "mem.hpp"

int main(void){

  double **** my4dptr = allocate4d(3,10,10,10);

  std::cout << my4dptr[2][3][4][5] << std::endl;
  my4dptr[2][3][4][5] = 3.0;
  std::cout << my4dptr[2][3][4][5] << std::endl;
}
