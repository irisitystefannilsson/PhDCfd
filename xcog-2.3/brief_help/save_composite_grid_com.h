#define ncom 3
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] ="ascii-format";
argument[0] = 1;
brief_help[0] = "Save the composite grid in an ASCII formatted file. This file "
"will be larger than the HDF-file, but it is readable by a (trained) human eye. "
"The format is described in the User's guide to Xcog, version 2.0, which is "
"available at http://www.na.chalmers.se/~andersp/xcog.html.";

command[1] ="hdf-format";
argument[1] = 1;
brief_help[1] = "Save the composite grid in a HDF file. This file will be "
"smaller than the ASCII file and take less time to read. The file is machine "
"independent and its contents can be examined with an HDF utility program "
"like hdp or vshow. The format is described in the User's guide to Xcog, "
"version 2.0, which is available at http://www.na.chalmers.se/~andersp/xcog.html.";

command[2] ="help";
argument[2] = 0;
brief_help[2] = NULL;

command[3] ="cancel";
argument[3] = 0;
brief_help[3] = "Don't save the overlapping grid on file and exit from this "
"command level.";

