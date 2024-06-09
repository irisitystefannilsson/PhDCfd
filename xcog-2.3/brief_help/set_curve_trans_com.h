#define ncom 6
  char *command[ncom+1], *brief_help[ncom+1];
  int argument[ncom+1], *save_on_copy=NULL;

command[0] ="rotation-angle"; argument[0] = 1;
command[1] ="horizontal-translation";  argument[1] = 1;
command[2] ="vertical-translation";  argument[2] = 1;

command[3] ="scale-factor";   argument[3] = 1;
command[4] ="show";           argument[4] = 0;
command[5] ="help";           argument[5] = 0;

command[6] ="exit";           argument[6] = 0;

brief_help[0] = "Set the rotation angle. The rotation is done before the "
"translation and the scaling";

brief_help[1] = "Set the horizontal translation distance. The translation is "
"done after the rotation, but before the scaling.";

brief_help[2] = "Set the vertical translation distance. The translation is "
"done after the rotation, but before the scaling.";

brief_help[3] = "Set the scale factor. The scaling is done after the rotation "
"and the translation.";

brief_help[4] = "Print the roation angle, translation distance and scale factor.";

brief_help[5] = NULL;

brief_help[6] = "Stop changing the transformation parameters and exit from this "
"command level.";
