/* Program to find number of multigrid levels
 * Run as 
 *    mglvls 257
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){
   int l=1, n;

   if(argc == 1){
      printf("Specify a grid size\n");
      exit(0);
   }

   n = atoi(argv[1]);
   printf("  Level    Size\n");
   printf("%5d %8d\n", l, n);

   while( (n-2)%2 != 0 ){
      n = (n-3)/2 + 2;
      l = l + 1;
      printf("%5d %8d\n", l, n);
   }

   return 0;

}
