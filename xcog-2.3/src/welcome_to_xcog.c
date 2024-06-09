#include <stdio.h>
#include "version.h"
extern int printf(const char *format, ...);

void welcome_to_xcog( void ){
  printf(
"\n"
"Welcome to Xcog version "VERSION".\n"
"Try the `tutorial' command for an introduction to Xcog and overlapping grids.\n"
"\n"
"A user's guide to the program can be found at\n"
"http://www.na.chalmers.se/~andersp/xcog/xcog.html\n"
"\n"
"Xcog was written by Anders Petersson.\n"
"The X-windows interface was written by Olof Runborg.\n"
"Please direct bug reports to andersp@nada.kth.se.\n"
"\n");
}

