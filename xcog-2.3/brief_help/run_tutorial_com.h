enum {STEPS, C_INTERP, CYL_IN_SQ, PAPER_KNIFE, SHOCK_CYL, INT_FLOW, HYP_WING,
#ifdef THE_GAR
      FOIL_BOX,
#endif
      C_GRID, MEX_GULF, MISC, HELP, EXIT};

const int LAST_COM = EXIT, LEVEL=1;

#ifdef THE_GAR
  char *COMMAND[13], *BRIEF_HELP[13];
  int *ARGUMENT=NULL, *SAVE_ON_COPY=NULL;
#else
  char *COMMAND[12], *BRIEF_HELP[12];
  int *ARGUMENT=NULL, *SAVE_ON_COPY=NULL;
#endif

COMMAND[STEPS] = "overlapping-grid-steps";
BRIEF_HELP[STEPS] = "Present a brief outline of how Xcog is used to make an "
"overlapping grid.";

COMMAND[C_GRID] = "c-grid";
BRIEF_HELP[C_GRID] = "Show the C-grid capability of the hyperbolic grid generator.";

COMMAND[C_INTERP] = "command-interpreter";
BRIEF_HELP[C_INTERP] = "Show the basic features of the command interpreter.";

COMMAND[CYL_IN_SQ] = "cylinder-in-square";
BRIEF_HELP[CYL_IN_SQ] = "Create an overlapping grid around a cylinder in a square.";

COMMAND[PAPER_KNIFE] = "paper-knife";
BRIEF_HELP[PAPER_KNIFE] = "Generate a grid for modelling the paper coating process.";

COMMAND[SHOCK_CYL] = "shock-cylinder";
BRIEF_HELP[SHOCK_CYL] = "Make a grid for computing the supersonic flow around a "
"cylinder.";

COMMAND[INT_FLOW] = "internal-flow";
BRIEF_HELP[INT_FLOW] = "Discretize a combustion chamber with an outflow pipe. "
"This problem exemplifies the use of mixed boundaries.";

COMMAND[HYP_WING] = "hyperbolic-wing";
BRIEF_HELP[HYP_WING] = "Demonstrate the hyperbolic grid generator by making "
"a grid around a NACA-0012 wing profile with a refinement around a shock.";

#ifdef THE_GAR
COMMAND[FOIL_BOX] = "wing-in-box";
BRIEF_HELP[FOIL_BOX] = "Make a grid around a NACA-0012 wing profile.";
#endif

COMMAND[C_GRID] = "c-grid";
BRIEF_HELP[C_GRID] = "Show the C-grid capability of the hyperbolic grid generator.";

COMMAND[MEX_GULF] = "mexican-gulf";
BRIEF_HELP[MEX_GULF] = "Discretize the gulf of Mexico by reading in externally "
"defined grid points in plot3d ASCII format.";

COMMAND[MISC] = "miscellaneous-features";
BRIEF_HELP[MISC] = "Demonstrate the following features of Xcog:\n"
"o  Reading in a saved state.\n"
"o  How to make postscript copies of the contents of the different windows.\n"
"o  How to modify the graphical contents of the different windows.";

COMMAND[HELP] = "help";
BRIEF_HELP[HELP] = NULL;

COMMAND[EXIT] = "exit";
BRIEF_HELP[EXIT] = "Exit from the tutorial mode and return to the top command "
"level.";


