enum{RED, GREEN, BLUE, BROWN, PURPLE, BLACK, GREY};

const int LAST_COM = GREY, LEVEL=1;
char *COMMAND[7], *BRIEF_HELP[7];
int *ARGUMENT=NULL, *SAVE_ON_COPY=NULL;

COMMAND[RED] = "red";
BRIEF_HELP[RED] = "Set the color to red.";

COMMAND[GREEN] = "green";
BRIEF_HELP[GREEN] = "Set the color to green.";

COMMAND[BLUE] = "blue";
BRIEF_HELP[BLUE] = "Set the color to blue.";

COMMAND[BROWN] = "brown";
BRIEF_HELP[BROWN] = "Set the color to brown.";

COMMAND[PURPLE] = "purple";
BRIEF_HELP[PURPLE] = "Set the color to purple.";

COMMAND[BLACK] = "black";
BRIEF_HELP[BLACK] = "Set the color to black.";

COMMAND[GREY] = "grey";
BRIEF_HELP[GREY] = "Set the color to grey.";

/* COMMAND[CANCEL] = "cancel"; */
/* BRIEF_HELP[CANCEL] = "Do not change the color"; */

