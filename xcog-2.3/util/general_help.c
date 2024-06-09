#include <stdio.h>
extern int printf(const char *format, ...);

void general_help( void ){
  printf(
"The user interface employs a tcsh-style command interpreter with the following features:\n"
" o  TAB completion of commands.\n"
" o  Upwards arrow or ^P traverses up in the command history list.\n"
" o  Downwards arrow or ^N goes down in the command history list.\n"
" o  Only the unique part between each `-' of a command needs to be given.\n"
" o  ? lists all possible commands.\n"
" o  !cmnd executes the shell command `cmnd'\n\n"
"The graphics windows have the following built in functionality:\n"
" o  Left mouse button zooms in.\n"
" o  Middle mouse button shows all.\n"
" o  Right mouse button zooms out.\n\n"
);
}
