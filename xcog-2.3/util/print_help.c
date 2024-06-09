#include <stdio.h>
extern int printf(const char *format, ...);

void print_help( char *brief_help ){
  int i, pos, length;
  const int min_length=70;

  if (brief_help == NULL) return;

  length = strlen(brief_help);

  pos = 0;
  printf("\n");
  do{
/* print one line that is at least 60 characters wide (if the string is
   sufficiently wide) and thereafter breaks the line at the first space.
*/
    for (i=0; i<min_length && pos < length && brief_help[pos] != '\n'; i++ )
      printf("%c", brief_help[pos++]);
/* wait with the newline until there is a space. Don't print the space */
    while (pos < length && brief_help[pos++] != ' ' && brief_help[pos-1] != '\n')
      printf("%c", brief_help[pos-1]);
    printf("\n");
  }
  while (pos<length);
}
	
