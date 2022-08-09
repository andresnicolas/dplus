
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct global_parameters {
   int    UseTab_wEoS;
   char   OutputFile[500];
   char   wEoSFile[500];
   double OmegaMatter;
   double OmegaDarkEnergy;
   double Hubble;
   double w0EoS;
   double waEoS;
   double InitialRedshift;
   double FinalRedshift;
   int    Nbins;   
};
extern struct global_parameters P;

void read_inputfile(char *filename);

