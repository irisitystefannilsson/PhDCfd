#include "xplot.h"
#include "real.h"

#define N_COLS 6

static color_type cc[] = {
   {0.0, 1.0, 0.0},
   {0.2, 0.0, 1.0},
   {0.9, 0.9, 0.0},
   {0.5, 0.0, 0.8},
   {1.0, 0.1, 0.0},
   {0.7, 0.6, 0.5}
 };

void main(){
  FILE *fp;
  input_output io, *io_ptr=&io;
  char *titles[] = {"Composite Grid", "Curve", "Grid"};
  real xy[][2] = {{0.0, 1.0}, {0.0, 0.35}, {0.35, 0.35}};
  real width[]  = {0.7, 0.34, 0.35};
  real height[] = {0.6, 0.3, 0.3};
  real xyab[2][2]= {-1.0, 1.0, -1.0, 1.0}, x;

  int i, font_number, h_center, v_center, w=0, ret;
  const int ncom=8;
  char *command[ncom];


  if (!PL_init(io_ptr, titles, xy, width, height, 3))
    printf("Warning! No X-windows interface.\nIt's OK, but you won't see"
	   " anything.\n");

  PL_scale(0, TRUE);
  PL_scale(1, TRUE);
/*  PL_window(0, xyab);
  PL_erase(0);
  PL_start_plot(0);
*/
  for (w=0; w<3; w++) {
/*    PL_window(w, xyab);*/
    PL_erase(w);
    PL_start_plot(w);
    PL_color(cc[w]);

    /* draw a cross */
    PL_move(-1.0,-1.0*w);
    PL_draw( 1.0, 1.0);
    PL_move(-1.0, 1.0);
    PL_draw( 1.0,-1.0);

    /* draw a number of boxes */
    for (i=1; i<=5; i++){
      PL_color(cc[(i*7) % N_COLS]);
      x = i/5.0;
      PL_move(-x,-x);
      PL_draw( x*w,-x);
      PL_draw( x, x);
      PL_draw(-x, x);
      PL_draw(-x,-x);
    }

    /* plot a couple of right adjusted strings */
    v_center = 1;
    h_center = 1-w;
    font_number = 3;
    PL_color(cc[(w*5) % N_COLS]);
    PL_plot_string("Hello world", 0.8, 0.7, h_center, v_center, font_number);
    PL_plot_string("Hello world", 0.8, 0.5, h_center, v_center, font_number);
    PL_plot_string("Hello world", 0.8, 0.3, h_center, v_center, font_number);
    PL_plot_string("Hello world", 0.8, 0.1, h_center, v_center, font_number);

    /* plot a couple of left adjusted strings */
    v_center = 1;
    h_center = -1;
    PL_plot_string("Hello world", 0.8,-0.1, h_center, v_center, font_number);
    PL_plot_string("Hello world", 0.8,-0.3, h_center, v_center, font_number);
    PL_plot_string("Hello world", 0.8,-0.5, h_center, v_center, font_number);
    PL_plot_string("Hello world", 0.8,-0.7, h_center, v_center, font_number);

    /* plot a couple of centered strings */
    v_center = 1-w;
    h_center = 0;
    PL_color(cc[(w*2) % N_COLS]);
    PL_plot_string("Hello world",-0.8, 0.7, h_center, v_center, font_number);
    PL_plot_string("Hello world",-0.8, 0.5, h_center, v_center, font_number);
    PL_plot_string("Hello world",-0.8, 0.3, h_center, v_center, font_number);
    PL_plot_string("Hello world",-0.8, 0.1, h_center, v_center, font_number);
    PL_plot_string("Hello world",-0.8,-0.1, h_center, v_center, font_number);
    PL_plot_string("Hello world",-0.8,-0.3, h_center, v_center, font_number);
    PL_plot_string("Hello world",-0.8,-0.5, h_center, v_center, font_number);
    PL_plot_string("Hello world",-0.8,-0.7, h_center, v_center, font_number);

    /* check the vertical centering */
    v_center = 1;
    h_center = 1;
    PL_plot_string("Hello world",0.0, 0.8, h_center, v_center, font_number);
    v_center = -1;
    h_center = -1;
    PL_plot_string("Hello world",0.0, 0.8, h_center, v_center, font_number);
    v_center = 1;
    h_center = 1;
    PL_plot_string("Hello world",0.0, 0.6, h_center, v_center, font_number);
    v_center = -1;
    h_center = -1;
    PL_plot_string("Hello world",0.0, 0.6, h_center, v_center, font_number);
    v_center = 1;
    h_center = w-1;
    PL_plot_string("Hello world",0.0, 0.4, h_center, v_center, font_number);
    v_center = -1;
    h_center = -1;
    PL_plot_string("Hello world",0.0, 0.4, h_center, v_center, font_number);
    v_center = 1;
    h_center = 1;
    PL_plot_string("Hello world",0.0, 0.2, h_center, v_center, font_number);
    v_center = -1;
    h_center = -1;
    PL_plot_string("Hello world",0.0, 0.2, h_center, v_center, font_number);

    /* check the horizontal placing of individual characters */
    v_center = 0; h_center = 0;

    PL_color(cc[(w*2) % N_COLS]);
    PL_marker(-0.8, -0.82, DOT+w); 
    PL_marker(-0.6, -0.6*w, CIRCLE+w);
    PL_marker(-0.4, -0.4, BOX+w);
    PL_color(cc[(w*15) % N_COLS]);
    PL_marker(-0.2, -0.2, TRIANGLE);
    PL_marker(-0.2,  0.2, FILLED_TRIANGLE);
    PL_marker( 0.2,  0.2, INV_TRIANGLE);
    PL_marker( 0.0,  0.0, FILLED_INV_TRIANGLE);
    PL_marker(0.4, -0.4, FILLED_CIRCLE+w);
    PL_marker(0.2, -0.2, FILLED_BOX+w);
    PL_marker(0.2, -0.4, CROSS);

  }
  PL_end_plot();

  command[0] = "exit";
  command[1] = "scale-off";
  command[2] = "scale-on";
  command[3] = "sesam-e-kul";
  command[4] = "title-on";
  command[5] = "title-off";
  command[6] = "dumpa-ps";
  command[7] = "cancel";

  PL_title(1, "Och s} sa han: kom nu ketchup s} g}r vi");
  PL_title(2, "Bakom stenen stod Urban och slet");

  while ((ret=get_command( io_ptr, "PLtest>", command, ncom-1, 0, NULL, NULL)) != 0 ) {
    if (ret == 1) PL_scale(0, FALSE);
    if (ret == 2) PL_scale(0, TRUE);
    if (ret == 4) PL_title(0, "Voffeltoffel");
    if (ret == 5) PL_title(0, NULL);
    if (ret == 6) {
      fp = fopen("wowow.ps", "w");
      PL_postscript(0, fp, 0);
      fclose(fp);
    }
  }

/* remove the windows */
  PL_reset( io_ptr );
}
