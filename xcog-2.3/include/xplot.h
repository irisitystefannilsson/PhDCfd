#ifndef xplot_h
#define xplot_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <termio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "real.h"
#include "input_output.h"

/* Defines for X inter routines */

#define MAX_VERTICES               100
#define MAX_WINDOWS                20
#define MAX_MOUSE_BUTTONS          5
#define MAX_POPUP_ENTRIES          20
#define RESET                      2
#define TRUE                       1
#define FALSE                      0
#define OK                         TRUE
#define ERROR                      FALSE
#define BORDER_WIDTH               3
#define LINE                       0
#define POLYGON                    1
#define GENERIC                    2
#define POPUP_ENTRY_HEIGHT         20
#define TEXT_SAVE_BUF_SIZE         50
#define LINE_SAVE_BUF_SIZE         500
#define MARKER_SAVE_BUF_SIZE       100
#define COLOR_SAVE_BUF_SIZE        16

#define ConvertCoordX(x,xmin,xmax,pixels) (((pixels)*((x)-(xmin))/((xmax)-(xmin)))+0.5)
#define ConvertCoordY(y,ymin,ymax,pixels) (((pixels)*((ymax)-(y))/((ymax)-(ymin)))+0.5)

/* Types for X inter */
					     
typedef struct _textsave {
  int hjust, vjust;
  real x, y;
  char *s;
  unsigned int cind;
} TextSaveEntry;

typedef struct _linesave {
  real x0, y0, x1, y1;
  unsigned int cind, lw, dash;
} LineSaveEntry;

typedef struct _markersave {
  real x, y;
  int type;
  unsigned int cind;
} MarkerSaveEntry;

typedef Window winid_type;

typedef struct Windata {
  winid_type window;
  GC gc;
  int col_index, line_width, dash;
  int width, height;
  int xpadl, xpadr, ypadt, ypadb;
  real xmin, xmax, ymin, ymax;
  real sxmin, sxmax, symin, symax;
  real aspect;
  real cx, cy;
  TextSaveEntry *ts;
  int ts_count;
  int ts_size;
  LineSaveEntry *ls;
  int ls_count;
  int ls_size;
  MarkerSaveEntry *ms;
  int ms_count;
  int ms_size;
} win_data_type;

typedef struct Entry {
  char *name;
  int ret_val;
  /*	  popup_type *next_popup;*/
} popup_entry;

typedef struct Popup {
  int n_entries;
  popup_entry *entry;
} popup_type;


/* Function prototypes for PL functions */

void
PL_start_dash(void);
void
PL_stop_dash(void);
void 
PL_set_same_scale( int window_number, int flag );
void 
PL_set_labels( int window_number, char *new_x_label, char *new_y_label );
void 
PL_scale(int, int);
void 
PL_title(int, char *);
int 
PL_init(input_output *io_ptr, char **title, real xy[][2], real *width, 
	real *height, int n_windows, int use_graphics);
void 
PL_reset(input_output *io_ptr);
void
PL_start_plot(int);
void 
PL_end_plot(void);
void 
PL_erase(int);
void 
PL_window(int, real[2][2]);
void 
PL_color(color_type);
void 
PL_move(real,real);
void 
PL_draw(real,real);
void 
PL_poll(void);
void 
PL_marker(real x, real y, int type);
void 
PL_plot_string(char *text, real x, real y, int h_center, int v_center, int font_size);
void 
PL_postscript(int w, FILE *fp, int col_flag);
void 
get_new_color( color_type *color_ptr, char type );
void
PL_line_width(unsigned int line_width);

/* Defines for command interpreter */

#define RAWON     1
#define RAWOFF    0
#define RAWNDELON 2

/* Command buffer sizes */

#define MAX_LINE_SIZE 256
#define SCROLLBACK_SIZE 100

/* Function prototypes for command interpreter */

int get_command(input_output *io_ptr, char *prompt,
		char **command, int n_command,
		int command_level, int *save_on_copy, int *argument);
char *get_word(input_output *io_ptr, char *prompt, char *deflt,
	       int save_on_copy);
int get_int(input_output *io_ptr, char *prompt, int deflt, int level);
real get_real(input_output *io_ptr, char *prompt, real deflt, int level);
void getline_no_block(char *prompt, char *buf, char **cmds, int n);
void wait_for_key(char *prompt);
void general_help( void );
void print_help( char *brief_help );
FILE *open_binary_file(input_output *io_ptr, char *prompt, char *deflt, 
		     char read_write, int file_type);
FILE *open_ascii_file( input_output *io_ptr, char *prompt, char *deflt, 
		      char **file_name, char read_write, int file_type, 
		      int save_on_copy);
FILE *open_this_ascii_file( char *file_name, char read_write, int file_type, 
			   int quiet);
FILE *open_this_binary_file( char *file_name, int file_type, int quiet);
int get_on_off( input_output *io_ptr, char *prompt );
int get_yes_no( input_output *io_ptr, char *prompt, int level );


#endif
