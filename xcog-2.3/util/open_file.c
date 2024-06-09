#include "xplot.h"
/*#include "stupid_compiler.h"*/

static void wrong_file_type( char info16 ){
  printf("Error: wrong type of file!\n");
  switch( info16 ){

  case '1':
    printf("This is an ascii command file\n");
    break;

  case '2':
    printf("This is a binary curve + mapping file\n");
    break;

  case '3':
    printf("This is a binary overlapping grid file\n");
    break;

  case '4':
    printf("This is an ascii node-point file\n");
    break;

  case '5':
    printf("This is an ascii overlapping grid file\n");
    break;

  default:
    ;

  }
}


FILE *open_binary_file(input_output *io_ptr, char *prompt, char *deflt, 
		     char read_write, int file_type){
  char *token, info[19], template[19], *xcog_home, file_name_2[250];
  FILE *fp;

  token = get_word( io_ptr, prompt, deflt, 1);

  sprintf(template, "#xcog-file-type-%i\n", file_type);

/* if the env variable XCOG_HOME is not set, look in the current directory */
  if (!(xcog_home = getenv("XCOG_HOME")))
    sprintf(file_name_2, "./%s", token);
  else
    sprintf(file_name_2, "%s/data/%s", xcog_home, token);

  if (read_write == 'w') {
    if ((fp = fopen( token, "wb" )) == NULL &&
	(fp = fopen( file_name_2, "wb" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", token);
    else if (file_type != 0){
/* write an identifier line at the top line */
      write( fileno(fp), template, 18*sizeof(char) );
    }
  }
  else if (read_write == 'r') {
    if ((fp = fopen( token, "rb" )) == NULL &&
	(fp = fopen( file_name_2, "rb" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", token);
    else if (file_type != 0){
      read( fileno(fp), info, 18*sizeof(char) );
      info[18]='\0';
/* error message */
      if (strcmp(template, info)!=0){
	wrong_file_type( info[16] );
	fclose( fp );
	fp = NULL;
      }
    }
  }
  else{
    printf("Unknown read-write status: %c\n", read_write);
    fp = NULL;
  }

  return fp;
}



FILE *
open_ascii_file( input_output *io_ptr, char *prompt, char *deflt, 
		char **file_name, char read_write, int file_type, 
		int save_on_copy){
  char info[19], template[19], *xcog_home, file_name_2[250];
  FILE *fp=NULL;
  
  *file_name = get_word( io_ptr, prompt, deflt, save_on_copy );
  
  sprintf(template, "#xcog-file-type-%i\n", file_type);

/* if the env variable XCOG_HOME is not set, look in the current directory */
  if (!(xcog_home = getenv("XCOG_HOME")))
    sprintf(file_name_2, "./%s", *file_name);
  else
    sprintf(file_name_2, "%s/data/%s", xcog_home, *file_name);

  if (read_write == 'w'){
    if ((fp = fopen( *file_name, "w" )) == NULL &&
	(fp = fopen( file_name_2, "w" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", *file_name);
    else if (file_type != 0){
/* write an identifier line at the top line */
      write( fileno(fp), template, 18*sizeof(char) );
    }
  }
  else if (read_write == 'r'){
    if ((fp = fopen( *file_name, "r" )) == NULL &&
	(fp = fopen( file_name_2, "r" )) == NULL)
      printf("Warning: Failed to open the file `%s'\n", *file_name);
    else if (file_type != 0){
      read( fileno(fp), info, 18*sizeof(char) );
      info[18]='\0';
      if (strcmp(template, info)!=0){
/* error message */
	wrong_file_type( info[16] );
	fclose( fp );
	fp = NULL;
      }
    }
  }
  return fp;
}


FILE *
open_this_ascii_file( char *file_name, char read_write, int file_type, int quiet){
  char info[19], template[19], *xcog_home, file_name_2[250];
  FILE *fp=NULL;
  
  sprintf(template, "#xcog-file-type-%i\n", file_type);

/* if the env variable XCOG_HOME is not set, look in the current directory */
  if (!(xcog_home = getenv("XCOG_HOME")))
    sprintf(file_name_2, "./%s", file_name);
  else
    sprintf(file_name_2, "%s/data/%s", xcog_home, file_name);

  if (read_write == 'w'){
    if ((fp = fopen( file_name, "w" )) == NULL &&
	(fp = fopen( file_name_2, "w" )) == NULL){
      if (!quiet) printf("Warning: Failed to open the file `%s'\n", file_name);
    }
    else if (file_type != 0){
/* write an identifier line at the top line */
      write( fileno(fp), template, 18*sizeof(char) );
    }
  }
  else if (read_write == 'r'){
    if ((fp = fopen( file_name, "r" )) == NULL &&
	(fp = fopen( file_name_2, "r" )) == NULL){
      if (!quiet) printf("Warning: Failed to open the file `%s'\n", file_name);
    }
    else if (file_type != 0){
      read( fileno(fp), info, 18*sizeof(char) );
      info[18]='\0';
      if (strcmp(template, info)!=0){
/* error message */
	if (!quiet) wrong_file_type( info[16] );
	fclose( fp );
	fp = NULL;
      }
    }
  }
  return fp;
}


FILE *open_this_binary_file( char *file_name, int file_type, int quiet){
  char info[19], template[19], *xcog_home, file_name_2[250];
  FILE *fp;
  
  sprintf(template, "#xcog-file-type-%i\n", file_type);

/* if the env variable XCOG_HOME is not set, look in the current directory */
  if (!(xcog_home = getenv("XCOG_HOME")))
    sprintf(file_name_2, "./%s", file_name);
  else
    sprintf(file_name_2, "%s/data/%s", xcog_home, file_name);

  if ((fp = fopen( file_name, "rb" )) == NULL &&
      (fp = fopen( file_name_2, "rb" )) == NULL){
    if (!quiet) printf("Warning: Failed to open the file `%s'\n", file_name);
  }
  else if (file_type != 0){
    read( fileno(fp), info, 18*sizeof(char) );
    info[18]='\0';
    if (strcmp(template, info)!=0){
/* error message */
      if (!quiet) wrong_file_type( info[16] );
      fclose( fp );
      fp = NULL;
    }
  }
  return fp;
}
