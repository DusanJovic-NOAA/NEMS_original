#include <stdio.h>
#include <sys/resource.h>

/* void print_memory(void) */
void print_memory(char *c)
{
    long Kbytes;
    double Mbytes;
    struct rusage RU;

    getrusage(RUSAGE_SELF, &RU);
 
    Kbytes = RU.ru_maxrss;

    Mbytes = ((double) Kbytes) / 1024.0;

/*    printf("Memory used = %.2f MBytes\n", Mbytes); */
/*    printf("Memory used = %.2f MBytes c = %s \n", Mbytes, c); */
    fprintf(stderr,"Memory used = %.4f MBytes c = %s \n", Mbytes, c);
  
    return;
}
