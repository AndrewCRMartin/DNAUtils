#include <stdio.h>
#include "bioplib/seq.h"

int main(int argc, char **argv);
int PatchAAs(char *aas, int numaa, char newaa);

int main(int argc, char **argv)
{
   char bases[] = 
      {
         'A','T','C','G'
      }
   ;
   char codon[4];
   int i,j,k,l,numaa,base;
   char aas[80];
   char newaa;
   int  total = 0, numcodons=0;
   
   codon[3] = '\0';
   /* Construct each possible codon */
   for(i=0; i<4; i++)
   {
      codon[0] = bases[i];
      for(j=0; j<4; j++)
      {
         codon[1] = bases[j];
         for(k=0; k<4; k++)
         {
            codon[2] = bases[k];

            /* Store the native amino acid */
            aas[0] = DNAtoAA(codon);
            if(aas[0] != 'X')
            {
               numcodons++;
               numaa = 1;
               
               /* Now change each base in turn to each other one */
               for(base=0; base<3; base++)
               {
                  char orig = codon[base];
                  for(l=0; l<4; l++)
                  {
                     codon[base] = bases[l];
                     newaa = DNAtoAA(codon);
                     numaa = PatchAAs(aas, numaa, newaa);
                  }
                  codon[base] = orig;
               }
               
               /* Print results */
               printf("Codon: %s; Native: %c; Alternates: %d\n",
                      codon, aas[0], numaa-1);
               total += numaa;
            }
         }
      }
   }

   printf("Average number of aas encoded by each codon and its variants: %.1f\n",
          (float)total/(float)numcodons);
   
   return(0);
}

int PatchAAs(char *aas, int numaa, char newaa)
{
   int i;
   if(newaa == 'X')
      return(numaa);

   for(i=0; i<numaa; i++)
   {
      if(aas[i] == newaa)
      {
         return(numaa);
      }
   }
   aas[numaa++] = newaa;
   return(numaa);
}
