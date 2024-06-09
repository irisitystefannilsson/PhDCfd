/****************************************************************************
  xplot.c

  Some X convenience functions.

  by Olof Runborg, NADA KTH, 91-02-04

  PL connection, redraw, zoom, etc. by Olof Runborg, UCLA Math, 95-02-10

  Command interpreter added by Anders Petersson, Math. Dept. UCLA, 95-02-07
****************************************************************************/

#include "real.h"
#include "xplot.h"
#include "min_max.h"
/*#include "stupid_compiler.h"*/

static Display *disp;
static Screen *screen;
static int screen_num;
static Colormap cmap;
static Window root_win;
static unsigned long black, white;

static int active_type;
static int active_pos;
static win_data_type windata[MAX_WINDOWS];
static int n_windows=0;
static XFontStruct *font_info;
static int X_is_running = FALSE;
static int X_Color_Server = FALSE;
static XRectangle xrect;
static int *same_scale;
static char **x_label, **y_label;

static XColor *cm;      /* 'My' colormap */
static int cm_count;
static int cm_size;

#ifndef M_LN10
#define	M_LN10		2.30258509299404568402
#endif

static real X_exp10( real x ){
  return exp( M_LN10 * x );
}

static int X_nint( real x ){
  return (int) (x + 0.5);
}

static int
conv_coord_x(win_data_type *w, real x)
{
  int ret, ww;

  ww = w->width - w->xpadl - w->xpadr-1;

  ret = w->xpadl + X_nint(ww *(x - w->xmin)/(w->xmax - w->xmin));

  return ret;
}

static int
conv_coord_y(win_data_type *w, real y)
{
  int ret, hh;

  hh = w->height - w->ypadt - w->ypadb-1;

  ret = w->ypadt + X_nint(hh *(w->ymax - y)/(w->ymax - w->ymin));

  return ret;
}

static real
conv_coord_x_inv(win_data_type *w, int x)
{
  int ww;
  real ret;

  ww = w->width - w->xpadl - w->xpadr;


  ret = w->xmin + (w->xmax - w->xmin) *(x - w->xpadl) / ww ;

  return ret;
}


static real
conv_coord_y_inv(win_data_type *w, int y)
{
  int hh;
  real ret;

  hh = w->height - w->ypadt - w->ypadb;

  ret = w->ymin + (w->ymax - w->ymin) * (w->height - w->ypadb - y) / hh;

  return ret;
}


static int
winid_exists(winid_type winid, int *pos){
  *pos=0;
  
  if (winid == (Window) NULL)
    return FALSE;
  if (n_windows > 0) {
    while ((*pos < n_windows) && (winid != windata[*pos].window))
      (*pos)++;
    
    return(*pos != n_windows);
  } else
    return(FALSE);
}


static int
color_exists(color_type color, int *cpos)
{
  XColor c;
  
  c.red = 65535.0 * color.r;
  c.green = 65535.0 * color.g;
  c.blue = 65535.0 * color.b;
 
  *cpos=0;
  
  if (cm_count > 0) {
    while ((*cpos < cm_count) && 
	   ((c.red  != cm[*cpos].red) ||
	    (c.green != cm[*cpos].green) ||
	    (c.blue  != cm[*cpos].blue)))
      (*cpos)++;
    
    return(*cpos != cm_count);
  } else
    return(FALSE);
}


static int
new_windata_pos(int *pos){
  *pos=0;

  if (n_windows == MAX_WINDOWS) {
    while ((*pos < n_windows) && (windata[*pos].window != (Window) NULL))
      (void) *pos++;
    return(*pos != n_windows);
  } else 
    *pos = n_windows++;

  return(OK);
}

static void
XDrawMarker(Display *this_disp, Drawable d, GC gc, int x, int y, int type)
{
  XPoint p[3];

  switch (type) {
  case DOT:
    XDrawPoint(this_disp, d, gc, x, y);
    break;
  case CIRCLE:
    XDrawArc(this_disp, d, gc, x-3, y-3, 6, 6, 0, 360*64);
    break;
  case BOX:
    XDrawRectangle(this_disp, d, gc, x-3, y-3, 6, 6);
    break;
  case ASTERISK:
    XDrawLine(this_disp, d, gc, x-3, y, x+3, y);
    XDrawLine(this_disp, d, gc, x, y-3, x, y+3);
    XDrawLine(this_disp, d, gc, x-2, y-2, x+2, y+2);
    XDrawLine(this_disp, d, gc, x-2, y+2, x+2, y-2);
    break;
  case CROSS:
    XDrawLine(this_disp, d, gc, x-4, y-4, x+4, y+4);
    XDrawLine(this_disp, d, gc, x-4, y+4, x+4, y-4);
    break;
  case FILLED_CIRCLE:
    XDrawArc(this_disp, d, gc,  x-3, y-3, 6, 6, 0, 360*64);
    XFillArc(this_disp, d, gc,  x-3, y-3, 6, 6, 0, 360*64);
    break;
  case FILLED_BOX:
    XFillRectangle(this_disp, d, gc, x-3, y-3, 7, 7);
    break;
  case TRIANGLE:
    XDrawLine(this_disp, d, gc, x, y-5, x+4, y+2);
    XDrawLine(this_disp, d, gc, x+4, y+2, x-4, y+2);
    XDrawLine(this_disp, d, gc, x-4, y+2, x, y-5);
    break;
  case FILLED_TRIANGLE:
    p[0].x = x;     p[0].y = y-6;
    p[1].x = x+5;   p[1].y = y+3;
    p[2].x = x-5;   p[2].y = y+3;
    XFillPolygon(this_disp, d, gc, p, 3, Convex, CoordModeOrigin);
    break;
  case INV_TRIANGLE:
    XDrawLine(this_disp, d, gc, x, y+5, x+4, y-2);
    XDrawLine(this_disp, d, gc, x+4, y-2, x-4, y-2);
    XDrawLine(this_disp, d, gc, x-4, y-2, x, y+5);
    break;
  case FILLED_INV_TRIANGLE:
    p[0].x = x;     p[0].y = y+6;
    p[1].x = x+5;   p[1].y = y-3;
    p[2].x = x-5;   p[2].y = y-3;
    XFillPolygon(this_disp, d, gc, p, 3, Convex, CoordModeOrigin);
    break;
  default:
    fprintf(stderr, "warning: marker type %d could not be displayed\n", type);
  }
}



static void
add_color_entry(color_type col)
{
  XColor *tmp, c;
  
  if (cm_count == cm_size) {
    tmp = (XColor *) malloc(sizeof(XColor)*cm_size*2);
    memcpy(tmp, cm, sizeof(XColor)*cm_size);  
    free(cm);
    cm = tmp;

    cm_size *= 2;
/*    printf("Colors: %d\n", cm_size);*/
  }

  cm[cm_count].red   = 65535.0 * col.r;
  cm[cm_count].green = 65535.0 * col.g;
  cm[cm_count].blue  = 65535.0 * col.b;

  cm[cm_count].pixel = black;
  if (cm[cm_count].red   == 65535 &&
      cm[cm_count].green == 65535 &&
      cm[cm_count].blue  == 65535)
    cm[cm_count].pixel = white;
  
  if (X_Color_Server) {
    c = cm[cm_count];
    if (XAllocColor(disp, cmap, &c))
      cm[cm_count].pixel = c.pixel;
  }

/*  fprintf(stderr, "col %d: %d %d %d %ld\n",
	  cm_count, 
	  cm[cm_count].red,
	  cm[cm_count].green,
	  cm[cm_count].blue,
	  cm[cm_count].pixel);*/

  cm_count++;
}

static void
add_textsave_entry(int pos, real x, real y, 
		   int hjust, int vjust, char *s)
{
  TextSaveEntry *tmp;

  if (windata[pos].ts_count == windata[pos].ts_size) {
    tmp = (TextSaveEntry *) malloc(sizeof(TextSaveEntry)*
				   windata[pos].ts_size*2);
    memcpy(tmp, windata[pos].ts, sizeof(TextSaveEntry)*windata[pos].ts_size);
    free(windata[pos].ts);
    windata[pos].ts = tmp;

    windata[pos].ts_size *= 2;
/*    printf("Text, Win: %d, %d\n", windata[pos].ts_size, pos);*/
  }

  windata[pos].ts[windata[pos].ts_count].x     = x;
  windata[pos].ts[windata[pos].ts_count].y     = y;
  windata[pos].ts[windata[pos].ts_count].hjust = hjust;
  windata[pos].ts[windata[pos].ts_count].vjust = vjust;
  windata[pos].ts[windata[pos].ts_count].s = (char *) malloc(strlen(s)+1);
  strcpy(windata[pos].ts[windata[pos].ts_count].s, s);
  
  windata[pos].ts[windata[pos].ts_count].cind = windata[pos].col_index;

  windata[pos].ts_count++;
}

static void
replay_textsave(int pos)
{
  int i, xi, yi, l, w;
  int ci=-1;

  for (i=0; i<windata[pos].ts_count; i++) {
    if (X_Color_Server && ci != windata[pos].ts[i].cind) {
      XSetForeground(disp, windata[pos].gc, cm[windata[pos].ts[i].cind].pixel);
      ci = windata[pos].ts[i].cind;
    }

    xi = conv_coord_x(windata+pos, windata[pos].ts[i].x);
    yi = conv_coord_y(windata+pos, windata[pos].ts[i].y);
    l = strlen(windata[pos].ts[i].s);
    w = XTextWidth(font_info, windata[pos].ts[i].s, l);
    
    if (windata[pos].ts[i].hjust == 0)   /* Centered */
      xi = xi - w / 2 + 1;
    if (windata[pos].ts[i].hjust == -1)   /* Right justified */
      xi = xi - w + 1;
    
    if (windata[pos].ts[i].vjust == 1)   /* Bottom justified */
      yi = yi + font_info->ascent;
    if (windata[pos].ts[i].vjust == 0)   /* Centered */
      yi = yi + (font_info->ascent + font_info->descent) / 2 - 1;
    if (windata[pos].ts[i].vjust == -1)   /* Top justified */
      yi = yi - font_info->descent + 1;
    
    XDrawString(disp, windata[pos].window, windata[pos].gc,
		xi, yi, windata[pos].ts[i].s, l);
    
  }
  if (X_Color_Server)
    XSetForeground(disp, windata[pos].gc, cm[windata[pos].col_index].pixel);
  XFlush(disp);
}

static void
reset_textsave(int pos){
  int i;

  for (i=0; i<windata[pos].ts_count; i++)
    free(windata[pos].ts[i].s);
  windata[pos].ts_count = 0;
}


static void
add_linesave_entry(int pos, real x0, real y0, real x1, real y1)
{
  LineSaveEntry *tmp;

  if (windata[pos].ls_count == windata[pos].ls_size) {
    tmp = (LineSaveEntry *) malloc(sizeof(LineSaveEntry)*
				   windata[pos].ls_size*2);
    memcpy(tmp, windata[pos].ls, sizeof(LineSaveEntry)*windata[pos].ls_size);
    free(windata[pos].ls);
    windata[pos].ls = tmp;

    windata[pos].ls_size *= 2;
/*    printf("Lines, Win: %d, %d\n", windata[pos].ls_size, pos);*/
  }

  windata[pos].ls[windata[pos].ls_count].x0 = x0;
  windata[pos].ls[windata[pos].ls_count].y0 = y0;
  windata[pos].ls[windata[pos].ls_count].x1 = x1;
  windata[pos].ls[windata[pos].ls_count].y1 = y1;
  windata[pos].ls[windata[pos].ls_count].dash = windata[pos].dash;
  windata[pos].ls[windata[pos].ls_count].lw = windata[pos].line_width;
  windata[pos].ls[windata[pos].ls_count].cind = windata[pos].col_index;
  windata[pos].ls_count++;
}

static void
replay_linesave(int pos)
{
  int i, xi0, yi0, xi1, yi1;
  int ci=-1;
  int line_width = -1;
  
  for (i=0; i<windata[pos].ls_count; i++) {
      
    if (!((windata[pos].ls[i].x0 < windata[pos].xmin && 
	   windata[pos].ls[i].x1 < windata[pos].xmin) ||
	  (windata[pos].ls[i].x0 > windata[pos].xmax && 
	   windata[pos].ls[i].x1 > windata[pos].xmax) ||
	  (windata[pos].ls[i].y0 < windata[pos].ymin && 
	   windata[pos].ls[i].y1 < windata[pos].ymin) ||
	  (windata[pos].ls[i].y0 > windata[pos].ymax && 
	   windata[pos].ls[i].y1 > windata[pos].ymax)))
      {
	if (X_Color_Server && ci != windata[pos].ls[i].cind) {
	  XSetForeground(disp, windata[pos].gc, cm[windata[pos].ls[i].cind].pixel);
	  ci = windata[pos].ls[i].cind;
	}
	if (line_width != windata[pos].ls[i].lw) {
	  XSetLineAttributes(disp, windata[pos].gc, windata[pos].ls[i].lw, 
			     LineSolid, CapButt, JoinMiter);
	  line_width = windata[pos].ls[i].lw;
	}
	xi0 = conv_coord_x(windata+pos, windata[pos].ls[i].x0);
	yi0 = conv_coord_y(windata+pos, windata[pos].ls[i].y0);
	xi1 = conv_coord_x(windata+pos, windata[pos].ls[i].x1);
	yi1 = conv_coord_y(windata+pos, windata[pos].ls[i].y1);

	XDrawLine(disp, windata[pos].window, windata[pos].gc, 
		  xi0, yi0, xi1, yi1);
      }
  }
  if (X_Color_Server)
    XSetForeground(disp, windata[pos].gc, cm[windata[pos].col_index].pixel);
  XSetLineAttributes(disp, windata[pos].gc, windata[pos].line_width, 
		     LineSolid, CapButt, JoinMiter);
  XFlush(disp);
}



static void
reset_linesave(int pos)
{
  windata[pos].ls_count = 0;
}





static void
add_markersave_entry(int pos, real x, real y, int t)
{
  MarkerSaveEntry *tmp;

  if (windata[pos].ms_count == windata[pos].ms_size) {
    tmp = (MarkerSaveEntry *) malloc(sizeof(MarkerSaveEntry)*
				   windata[pos].ms_size*2);
    memcpy(tmp, windata[pos].ms, sizeof(MarkerSaveEntry)*windata[pos].ms_size);
    free(windata[pos].ms);
    windata[pos].ms = tmp;

    windata[pos].ms_size *= 2;
/*    printf("Markers, Win: %d, %d\n", windata[pos].ms_size, pos);*/
  }

  windata[pos].ms[windata[pos].ms_count].x = x;
  windata[pos].ms[windata[pos].ms_count].y = y;
  windata[pos].ms[windata[pos].ms_count].type = t;
  windata[pos].ms[windata[pos].ms_count].cind = windata[pos].col_index;
  windata[pos].ms_count++;
}

static void
replay_markersave(int pos)
{
  int i, x, y;
  int ci=-1;

  for (i=0; i<windata[pos].ms_count; i++) {
      
    if (!(windata[pos].ms[i].x < windata[pos].xmin ||
	  windata[pos].ms[i].x > windata[pos].xmax ||
	  windata[pos].ms[i].y < windata[pos].ymin ||
	  windata[pos].ms[i].y > windata[pos].ymax ) )
      {
	if (X_Color_Server && ci != windata[pos].ms[i].cind) {
	  XSetForeground(disp, windata[pos].gc, cm[windata[pos].ms[i].cind].pixel);
	  ci = windata[pos].ms[i].cind;
	}
	x = conv_coord_x(windata+pos, windata[pos].ms[i].x);
	y = conv_coord_y(windata+pos, windata[pos].ms[i].y);
	XDrawMarker(disp, windata[pos].window, windata[pos].gc, 
		    x, y, windata[pos].ms[i].type);
      }
  }
  if (X_Color_Server)
    XSetForeground(disp, windata[pos].gc, cm[windata[pos].col_index].pixel);
  XFlush(disp);
}



static void
reset_markersave(int pos)
{
  windata[pos].ms_count = 0;
}




static int
X_winaspect(winid_type winid, real *aspect){
  XWindowAttributes attrib;
  int pos;
  
  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {
    XGetWindowAttributes(disp, winid, &attrib);
    * aspect = ((real) (attrib.width-
			  windata[pos].xpadl-windata[pos].xpadr))
      / ((real) (attrib.height-
		   windata[pos].ypadt-windata[pos].ypadb));
    return(OK);
  } else
    return(ERROR);
}


static int X_poll_zoom(winid_type *winid, 
	    real *x0d, real *y0d, real *x1d, real *y1d){
  XEvent event;
  int pos, x0=0, y0=0, x1, y1, x1f, y1f, tmp, dummy;
  unsigned int mask;
  real x0t, y0t, x1t, y1t;

  if (!X_is_running)
    return ERROR;

  if (XCheckMaskEvent(disp, ButtonPressMask, &event)) {
    if (event.type == ButtonPress) {
      if (!XCheckMaskEvent(disp, ButtonReleaseMask, &event)) {
	*winid = event.xbutton.window;
	if (!winid_exists(*winid, &pos))
	  return ERROR;
	while (XCheckWindowEvent(disp, *winid, PointerMotionMask, &event)) {
	  if (event.type == MotionNotify) { 
	    x0 = event.xmotion.x;
	    y0 = event.xmotion.y;
	  }
	}
	x1=x0; y1=y0;
	x1f=x0; y1f=y0;
	XQueryPointer(disp, *winid, (Window *) &dummy, (Window *) &dummy, 
                  (int *) &dummy, (int *) &dummy, &dummy, &dummy, &mask);

	if ((mask & (Button1Mask |Button2Mask | Button3Mask)) == 
	    Button2Mask) { /* Reset */
	  while (!XCheckMaskEvent(disp, ButtonReleaseMask, &event))
	    ;
	  return RESET;
	}
        pos++; /* bug fix which I do not understand: Stefan */
	XSetFunction(disp, windata[pos].gc, GXinvert); 
	while (!XCheckMaskEvent(disp, ButtonReleaseMask, &event)) {
	  while (XCheckWindowEvent(disp, *winid, PointerMotionMask, &event)) {
	    if (event.type == MotionNotify) { 
	      x1 = event.xmotion.x;
	      y1 = event.xmotion.y;
	    }
	  }
	  if (x1 != x1f || y1 != y1f) {

	    /* Draw rubberband */
	    XDrawLine(disp, *winid, windata[pos].gc, x0, y0, x0, y1f);
	    XDrawLine(disp, *winid, windata[pos].gc, x0, y1f, x1f, y1f);
	    XDrawLine(disp, *winid, windata[pos].gc, x1f, y1f, x1f, y0);
	    XDrawLine(disp, *winid, windata[pos].gc, x1f, y0, x0, y0);
	    XDrawLine(disp, *winid, windata[pos].gc, x0, y0, x0, y1);
	    XDrawLine(disp, *winid, windata[pos].gc, x0, y1, x1, y1);
	    XDrawLine(disp, *winid, windata[pos].gc, x1, y1, x1, y0);
	    XDrawLine(disp, *winid, windata[pos].gc, x1, y0, x0, y0);
	    XFlush(disp);
	  }
	  x1f=x1; y1f=y1;
	}

	/* Erase rubberband */

	XDrawLine(disp, *winid, windata[pos].gc, x0, y0, x0, y1);
	XDrawLine(disp, *winid, windata[pos].gc, x0, y1, x1, y1);
	XDrawLine(disp, *winid, windata[pos].gc, x1, y1, x1, y0);
	XDrawLine(disp, *winid, windata[pos].gc, x1, y0, x0, y0);

	XFlush(disp);
	XSetFunction(disp, windata[pos].gc, GXcopy);

	/* Minimum size of zoom rectangle = 5x5 */
	
	if (winid_exists(*winid, &pos) && abs(x1-x0)>5 && abs(y1-y0)>5) {
	  if (x1 < windata[pos].xpadl || x1 > windata[pos].width-windata[pos].xpadr ||
	      y1 < windata[pos].ypadt || y1 > windata[pos].height-windata[pos].ypadb)
	    return FALSE;
	    
	  if (x0 > x1) { tmp=x1; x1=x0; x0=tmp; }
	  if (y1 > y0) { tmp=y1; y1=y0; y0=tmp; }
	    
	  if ((mask & (Button1Mask |Button2Mask | Button3Mask)) == 
	      Button1Mask ) { /* Zoom in */
	    *x0d = conv_coord_x_inv(windata+pos, x0);
	    *x1d = conv_coord_x_inv(windata+pos, x1);
	    *y0d = conv_coord_y_inv(windata+pos, y0);
	    *y1d = conv_coord_y_inv(windata+pos, y1);
	  }
	  if ((mask & (Button1Mask |Button2Mask | Button3Mask)) == 
		Button3Mask) { /* Zoom out */
	    x0t = conv_coord_x_inv(windata+pos, x0);
	    x1t = conv_coord_x_inv(windata+pos, x1);
	    y0t = conv_coord_y_inv(windata+pos, y0);
	    y1t = conv_coord_y_inv(windata+pos, y1);
	    *x0d = (windata[pos].xmin*(x1t-windata[pos].xmin)-
		    windata[pos].xmax*(x0t-windata[pos].xmin))/(x1t-x0t);
	    *x1d = (windata[pos].xmin*(windata[pos].xmax-windata[pos].xmin)-
		    (windata[pos].xmax-x0t)*(*x0d))/(x0t-windata[pos].xmin);
	    *y0d = (windata[pos].ymin*(y1t-windata[pos].ymin)-
		    windata[pos].ymax*(y0t-windata[pos].ymin))/(y1t-y0t);
	    *y1d = (windata[pos].ymin*(windata[pos].ymax-windata[pos].ymin)-
		    (windata[pos].ymax-y0t)*(*y0d))/(y0t-windata[pos].ymin);
	  }
	  return TRUE;
	}
      }
    }
  }
  return FALSE;
}

static int X_poll_resize(winid_type *winid){
  XEvent event;
  int pos;

  if (!X_is_running)
    return ERROR;

  if (XCheckMaskEvent(disp, StructureNotifyMask, &event))
    if (event.type == ConfigureNotify) {
      *winid = event.xconfigure.window;
      if (winid_exists(*winid, &pos)) {
	if (windata[pos].width != event.xconfigure.width ||
	    windata[pos].height != event.xconfigure.height) {
	  windata[pos].width = event.xconfigure.width;
	  windata[pos].height = event.xconfigure.height;
	  windata[pos].aspect = 
	    (real) (event.xconfigure.width-windata[pos].xpadl-windata[pos].xpadr) / 
	    (real) (event.xconfigure.height-windata[pos].ypadt-windata[pos].ypadb);
	  return TRUE;
	}
      }
    }
  return FALSE;
}


static int 
X_activate(winid_type winid)
{
  int pos;

  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {
    active_pos   = pos;
    active_type  = GENERIC;
    return(OK);
  } else
    return(ERROR);
}


static int 
X_color(color_type color)
{
  int cpos;

  if (!X_is_running)
    return ERROR;

/*fprintf(stderr, "Provar: %f %f %f\n", color.r, color.g, color.b);  */
  if (!color_exists(color, &cpos)) {
    add_color_entry(color);
    cpos = cm_count-1;
  }

  if (X_Color_Server)
    XSetForeground(disp, windata[active_pos].gc, cm[cpos].pixel);

  windata[active_pos].col_index = cpos;
  return(OK);
}

    
static void
X_lineto(real x1, real y1)
{
  int xi1, yi1, xi2, yi2;

  if (!X_is_running)
    return;

  if (!((x1 < windata[active_pos].xmin && windata[active_pos].cx < windata[active_pos].xmin) ||
      (x1 > windata[active_pos].xmax && windata[active_pos].cx > windata[active_pos].xmax) ||
      (y1 < windata[active_pos].ymin && windata[active_pos].cy < windata[active_pos].ymin) ||
	(y1 > windata[active_pos].ymax && windata[active_pos].cy > windata[active_pos].ymax)))
    {
      
      xi1 = conv_coord_x(windata+active_pos, windata[active_pos].cx);
      yi1 = conv_coord_y(windata+active_pos, windata[active_pos].cy);
      xi2 = conv_coord_x(windata+active_pos, x1);
      yi2 = conv_coord_y(windata+active_pos, y1);

      XDrawLine(disp, windata[active_pos].window, windata[active_pos].gc, 
		xi1, yi1, xi2, yi2);
      XFlush(disp);
    }

  add_linesave_entry(active_pos, windata[active_pos].cx, windata[active_pos].cy, x1, y1);

  windata[active_pos].cx=x1; windata[active_pos].cy=y1;
}


static void
X_marker(int type)
{
  int x=0, y=0;

  if (!X_is_running)
    return;

  if (!((windata[active_pos].cx < windata[active_pos].xmin) ||
        (windata[active_pos].cx > windata[active_pos].xmax) ||
        (windata[active_pos].cy < windata[active_pos].ymin) ||
	(windata[active_pos].cy > windata[active_pos].ymax)))
    {
      
      x = conv_coord_x(windata+active_pos, windata[active_pos].cx);
      y = conv_coord_y(windata+active_pos, windata[active_pos].cy);

      XDrawMarker(disp, windata[active_pos].window, 
		  windata[active_pos].gc, x, y, type);
      XFlush(disp);
    }

  add_markersave_entry(active_pos, windata[active_pos].cx, 
		       windata[active_pos].cy, type);

  windata[active_pos].cx=x; windata[active_pos].cy=y; 
}

static void
X_moveto(real x0, real y0)
{
  if (!X_is_running)
    return;

  windata[active_pos].cx=x0; windata[active_pos].cy=y0;
}


static void
X_text(char *s, int hjust, int vjust)
{    
  real xi, yi;
  int w, l;
  
  if (!X_is_running)
    return;

  xi = conv_coord_x(windata+active_pos, windata[active_pos].cx);
  yi = conv_coord_y(windata+active_pos, windata[active_pos].cy);

  l = strlen(s);
  w = XTextWidth(font_info, s, l);

  if (hjust == 0)   /* Centered */
    xi = xi - w / 2 + 1;
  if (hjust == -1)   /* Right justified */
    xi = xi - w + 1;

  if (vjust == 1)   /* Bottom justified */
    yi = yi + font_info->ascent;
  if (vjust == 0)   /* Centered */
    yi = yi + (font_info->ascent + font_info->descent) / 2 - 1;
  if (vjust == -1)   /* Top justified */
    yi = yi - font_info->descent + 1;

  XDrawString(disp, windata[active_pos].window, windata[active_pos].gc, xi, yi, s, l);
  XFlush(disp);

  add_textsave_entry(active_pos, windata[active_pos].cx, windata[active_pos].cy,
		     hjust, vjust, s);
}


static int
X_clear(winid_type winid)
{
  int pos;

  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {
    XClearWindow(disp, winid);
    reset_textsave(pos);
    reset_linesave(pos);
    reset_markersave(pos);
    return(OK);
  } else
    return(ERROR);
}

static int 
X_rescale(winid_type winid, real xmin, real ymin, real xmax, real ymax, int flag,
	  int same_scale_in_x_and_y)
/* if flag == FALSE = don't care about xmin... etc */
{
  real aspect;
  const real min_size=1.e-7;
  XWindowAttributes attrib;
  int pos;

  if (!X_is_running)
    return ERROR;

  if ((xmin <= xmax && ymin <= ymax) || flag == FALSE) {
    if (winid_exists(winid, &pos)) {
      XGetWindowAttributes(disp, winid, &attrib);
      windata[pos].width  = attrib.width;
      windata[pos].height = attrib.height;

      aspect = ((real) (attrib.width-
			  windata[pos].xpadl-windata[pos].xpadr))
		/ ((real) (attrib.height-
			     windata[pos].ypadt-windata[pos].ypadb));

      windata[pos].aspect = aspect;


      if (flag == FALSE) {
	xmin = windata[pos].sxmin;
	ymin = windata[pos].symin;
	xmax = windata[pos].sxmax;
	ymax = windata[pos].symax;
      }

      windata[pos].sxmin = xmin;
      windata[pos].symin = ymin;
      windata[pos].sxmax = xmax;
      windata[pos].symax = ymax;
      
      if (same_scale_in_x_and_y){

	if (aspect > (xmax-xmin) / real_max(min_size, ymax-ymin)) {
	  xmin = (xmin+xmax)/2 - real_max(min_size, ymax-ymin)*aspect/2;
	  xmax = xmin + real_max(min_size, ymax-ymin)*aspect;
	} else {
	  ymin = (ymin+ymax)/2 - real_max(min_size, xmax-xmin)/aspect/2;
	  ymax = ymin + real_max(min_size, xmax-xmin)/aspect;
	}
      
      }

      windata[pos].xmin = xmin;
      windata[pos].ymin = ymin;
      windata[pos].xmax = xmax;
      windata[pos].ymax = ymax;
      

      return(OK);
    }
  }
  return(ERROR);
}


/* x, y, w, h between 0 and 1! */

static winid_type 
X_winopen(char *name, real x, real y, real w, real h)
{
  XSetWindowAttributes attrib;
  XSizeHints size_hints;
  XEvent event;
  winid_type  new_winid;
  unsigned long valuemask;
  real aspect;
  int pos;

  if (!X_is_running)
    return (Window) NULL;

  if (new_windata_pos(&pos) == OK) {
  
    
    if (X_winaspect(root_win, &aspect)) {

	size_hints.x = WidthOfScreen(screen) * x;
	size_hints.y = HeightOfScreen(screen) * (1-y);
	size_hints.height = HeightOfScreen(screen) * h;
	size_hints.width = WidthOfScreen(screen) * w;

	size_hints.min_width = WidthOfScreen(screen) / 4;
	size_hints.min_height = size_hints.min_width / aspect;
	size_hints.flags = USPosition|USSize|PMinSize;

	attrib.background_pixel = white;
	attrib.border_pixel     = black;
	attrib.backing_store    = Always;
	attrib.event_mask       = EnterWindowMask|LeaveWindowMask|
	  ExposureMask|ButtonPressMask|ButtonReleaseMask|PointerMotionMask|
	    OwnerGrabButtonMask|StructureNotifyMask;
	
	valuemask = CWBackPixel|CWBorderPixel|CWBackingStore|CWEventMask;
	
	new_winid = XCreateWindow(disp, 
				  root_win, 
				  size_hints.x, 
				  size_hints.y,
				  size_hints.width, 
				  size_hints.height,
				  BORDER_WIDTH, 
				  CopyFromParent,
				  CopyFromParent,
				  CopyFromParent,
				  valuemask,
				  &attrib);
	
	XSetWMNormalHints(disp, new_winid, &size_hints);
	
	XStoreName(disp, new_winid, name);
	
	XMapWindow(disp, new_winid);
	XFlush(disp);
	
	/* Wait until window appear on screen */
	
	XWindowEvent(disp, new_winid, ExposureMask, &event);
	
	windata[pos].window = new_winid;

	windata[pos].gc  = XCreateGC(disp, root_win, 0, NULL);
	windata[pos].col_index = 1; /* Black as defined in X_init */
	windata[pos].line_width = 1; /* Thin lines */
	windata[pos].dash = FALSE; /* Solid lines */
	windata[pos].xpadl = 0;
	windata[pos].xpadr = 0;
	windata[pos].ypadt = 0;
	windata[pos].ypadb = 0;
	XSetForeground(disp, windata[pos].gc, black);
	XSetBackground(disp, windata[pos].gc, white);
	XSetFont(disp, windata[pos].gc, font_info->fid);

	windata[pos].ts = (TextSaveEntry *) 
	  malloc(sizeof(TextSaveEntry)*TEXT_SAVE_BUF_SIZE);
	windata[pos].ts_count = 0;
	windata[pos].ts_size = TEXT_SAVE_BUF_SIZE;

	windata[pos].ls = (LineSaveEntry *) 
	  malloc(sizeof(LineSaveEntry)*LINE_SAVE_BUF_SIZE);
	windata[pos].ls_count = 0;
	windata[pos].ls_size = LINE_SAVE_BUF_SIZE;

	windata[pos].ms = (MarkerSaveEntry *) 
	  malloc(sizeof(MarkerSaveEntry)*MARKER_SAVE_BUF_SIZE);
	windata[pos].ms_count = 0;
	windata[pos].ms_size = MARKER_SAVE_BUF_SIZE;

	X_rescale(new_winid, -1.0, -1.0, 1.0, 1.0, TRUE, 1);

	return(new_winid);
      } 
  }
  return(ERROR);
}


static int
X_winclose(winid_type winid)
{
  int pos;
  
  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {
    XDestroyWindow(disp, winid);
    windata[pos].window = (Window) NULL;
    reset_textsave(pos);
    free(windata[pos].ts);
    reset_linesave(pos);
    free(windata[pos].ls);
    reset_markersave(pos);
    free(windata[pos].ms);
    return(OK);
  } else
    return(ERROR);
}


static int 
X_init( void )
{
  char *display_name = NULL;
  Visual *default_visual;
  int default_depth;
  XVisualInfo vis_info;
  color_type c;
  int i;

  X_is_running = FALSE;
  X_Color_Server = FALSE;

  if ((disp = XOpenDisplay(display_name)) != 0) {
    
    root_win = DefaultRootWindow(disp);
    n_windows = 1;
    windata[0].window = root_win; 
    screen   = XDefaultScreenOfDisplay(disp);
    black    = BlackPixelOfScreen(screen);
    white    = WhitePixelOfScreen(screen);
    if ((font_info = XLoadQueryFont(disp, "8x13")) == NULL)
      font_info = XQueryFont(disp, XGContextFromGC(DefaultGCOfScreen(screen)));

    default_depth  = DefaultDepth(disp, screen_num);
    default_visual = DefaultVisual(disp, screen_num);
    cmap           = DefaultColormap(disp, screen_num);
    
    i=5; 
    while (!XMatchVisualInfo(disp, screen_num, default_depth,
			     i--, &vis_info))
      ;
    if (i >= StaticColor) {
      X_Color_Server = TRUE;
/*      fprintf(stderr, "xcog: color server detected.\n");*/
    }

    X_is_running = TRUE;

    cm = (XColor *)  malloc(sizeof(XColor)*COLOR_SAVE_BUF_SIZE);
    cm_count = 0;
    cm_size = COLOR_SAVE_BUF_SIZE;
    c.r=1.0; c.g=1.0; c.b=1.0; add_color_entry(c);  /* White */
    c.r=0.0; c.g=0.0; c.b=0.0; add_color_entry(c);  /* Black */

    return(OK);
    
  } else
    
    return(ERROR);
}


static void
X_exit( void )
{
  int i;

  if (!X_is_running)
    return;
  
  for (i=1; i<MAX_WINDOWS; i++)
    X_winclose(windata[i].window);

  free(cm);

  XCloseDisplay(disp);
  X_is_running = FALSE;
}


static int
X_redraw(winid_type winid)
{
  int pos;

  if (!X_is_running)
    return ERROR;
  
  if (winid_exists(winid, &pos)) {
    XClearWindow(disp, winid);
    replay_textsave(pos);
    replay_linesave(pos);
    replay_markersave(pos);
    return(OK);
  } else
    return(ERROR);
}


static int
X_winborder(winid_type winid, int xpadl, int xpadr, int ypadt, int ypadb, int flag)
     /* flag = FALSE -> don't care about pads */
{
  int pos;

  if (!X_is_running)
    return ERROR;
  
  if (winid_exists(winid, &pos)) {
    if (flag == TRUE) {
      windata[pos].xpadl = xpadl;
      windata[pos].xpadr = xpadr;
      windata[pos].ypadt = ypadt;
      windata[pos].ypadb = ypadb;
    }
    if (windata[pos].xpadl == 0 && 
	windata[pos].xpadr == 0 && 
	windata[pos].ypadt == 0 && 
	windata[pos].ypadb == 0)
      XSetClipMask(disp, windata[pos].gc, None);
    else {
      xrect.x = 0;
      xrect.y = 0;
      xrect.width  = windata[pos].width  - 
	windata[pos].xpadl -windata[pos].xpadr; 
      xrect.height = windata[pos].height - 
	windata[pos].ypadt - windata[pos].ypadb;
      XSetClipRectangles(disp, windata[pos].gc, 
			 windata[pos].xpadl, 
			 windata[pos].ypadt, 
			 &xrect, 1, Unsorted);
    }
  } else
    return(ERROR);
  return OK;
}

/***************************************
 * X scale functions 
 **************************************/


static void
tick_dist(real xmin, real xmax, real *tick1, real *dtick, int *ndec)
{
  real p, dnorm;

/* Compute a reasonable tick distribution on the axes */
  
  *ndec = (int) ( ceil(-log10(xmax-xmin)) + 0.5 );
  p = X_exp10(-(real) *ndec);
  dnorm = (xmax-xmin) / p;

  if (dnorm <= 3) 
    *dtick = 0.2; 
  else if (dnorm <= 6) 
    *dtick = 0.5; 
  else
    *dtick = 1.0;

  if (*dtick <  0.99) (*ndec)++ ;
  *dtick *= p;
  *tick1 = (floor(xmin / *dtick)+1)* (*dtick);

  if ((*ndec) < 0)
    *ndec = 0;

}

static void
get_scale_pad(int pos, 
	      int *xpadl, int *xpadr, int *ypadt, int *ypadb)
{
  real ymin, ymax, tick1, dtick;
  int ndec, l;
  char s[80];

  *xpadr = 0;
  *ypadt = 0;
  *ypadb = 20;
  
  ymin = conv_coord_y_inv(windata+pos, windata[pos].height);
  ymax = conv_coord_y_inv(windata+pos, 0);
  
  tick_dist(ymin, ymax, &tick1, &dtick, &ndec);

  if (ndec <= 4 && fabs(tick1) < 1.e4)
    sprintf(s, "%-.*f", ndec, tick1);
  else
    sprintf(s, "%-.3e", tick1);

  l = strlen(s);
  *xpadl = XTextWidth(font_info, s, l)+8;

/* number of ticks */
/*   n = X_nint((ymax-tick1)/dtick);  */

/* also compute the size of the last label */
/*   sprintf(s, "%-.*f", ndec, tick1+(n-1)*dtick); */
/*   l = strlen(s); */
/*   if (XTextWidth(font_info, s, l)+8 >*xpadl) */
/*     *xpadl=XTextWidth(font_info, s, l)+8; */
}


static int
X_scale_active(winid_type winid, int yesno)
{
  int pos;
  int xpadl, xpadr, ypadt, ypadb;

  if (!X_is_running)
    return ERROR;
  if (winid_exists(winid, &pos)) {
    if (yesno == 0) 
      X_winborder(winid, 0, 0, 0, 0, TRUE);
    else {
      get_scale_pad(pos, &xpadl, &xpadr, &ypadt, &ypadb);
      X_winborder(winid, xpadl+5, xpadr, ypadt, ypadb+5, TRUE);
    }
    return OK;
  }
  return ERROR;
}


static int
X_draw_scale(winid_type winid) 
{
  int pos, l, ndec;
  int x1, y1;
  int bordx, bordy, dummy;
  real xmin, xmax, ymin, ymax, t;
  real dtick, tick1;
  char s[80];
  XGCValues gcv;
  
  if (!X_is_running)
    return ERROR;
  
  if (winid_exists(winid, &pos)) {

    if (X_Color_Server) {
      XGetGCValues(disp, windata[pos].gc, GCForeground, &gcv);
      XSetForeground(disp, windata[pos].gc, black);
    }

/* No clip mask while drawing scales */

    XSetClipMask(disp, windata[pos].gc, None);

/* compute the size of the coordinate labels */
    get_scale_pad(pos, &bordx, &dummy, &dummy, &bordy);

/* Coordinate axis */

    XDrawLine(disp, windata[pos].window, windata[pos].gc, 
	      bordx,
	      windata[pos].height-bordy,
	      windata[pos].width-bordx,
	      windata[pos].height-bordy);
    XDrawLine(disp, windata[pos].window, windata[pos].gc, 
	      windata[pos].width-bordx,
	      windata[pos].height-bordy,
	      windata[pos].width-bordx-3,
	      windata[pos].height-bordy-3);
    XDrawLine(disp, windata[pos].window, windata[pos].gc, 
	      windata[pos].width-bordx,
	      windata[pos].height-bordy,
	      windata[pos].width-bordx-3,
	      windata[pos].height-bordy+3);

    l = XTextWidth(font_info, x_label[pos-1], strlen(x_label[pos-1]));
/*    x1 = windata[pos].width-bordx+1;*/
/*    y1 = windata[pos].height-bordy + (font_info->ascent + font_info->descent) / 2 - 1;*/
    x1 = windata[pos].width-bordx-l;
    y1 = windata[pos].height-bordy - (font_info->ascent + font_info->descent);
    XDrawString(disp, windata[pos].window, windata[pos].gc, x1, y1, x_label[pos-1], 
		strlen(x_label[pos-1]));

    XDrawLine(disp, windata[pos].window, windata[pos].gc, 
	      bordx, windata[pos].height-bordy,
	      bordx, bordy);
    XDrawLine(disp, windata[pos].window, windata[pos].gc, 
	      bordx, bordy, bordx-3, bordy+3);
    XDrawLine(disp, windata[pos].window, windata[pos].gc, 
	      bordx, bordy, bordx+3, bordy+3);

    l = XTextWidth(font_info, y_label[pos-1], strlen(y_label[pos-1]));
    x1 = bordx-l/2+1;
    y1 = bordy - 0.5*(font_info->ascent + font_info->descent);
    XDrawString(disp, windata[pos].window, windata[pos].gc, x1, y1, y_label[pos-1], 
		strlen(y_label[pos-1]));

/* Y-axis tick marks and labels */

    ymin = conv_coord_y_inv(windata+pos, windata[pos].height-bordy);
    ymax = conv_coord_y_inv(windata+pos, bordy);

    tick_dist(ymin, ymax, &tick1, &dtick, &ndec);

    for (t=tick1; t<ymax; t+=dtick) {
      y1 = conv_coord_y(windata+pos, t);
      XDrawLine(disp, windata[pos].window, windata[pos].gc,  /* Tick mark */
		bordx+3, y1, bordx-3, y1);
      if (fabs(t)< fabs(dtick)/2)
	strcpy(s, "0");
      else{
	if (ndec <= 4 && fabs(t) < 1.e4)
	  sprintf(s, "%-.*f", ndec, t);
	else
	  sprintf(s, "%-.3e", t);
      }
      l=strlen(s);
      y1 = y1 + (font_info->ascent + font_info->descent) / 2 - 1;
      XDrawString(disp, windata[pos].window, windata[pos].gc, 2, y1, s, l);
    }
    
    xmin = conv_coord_x_inv(windata+pos, bordx);
    xmax = conv_coord_x_inv(windata+pos, windata[pos].width-bordx);

    tick_dist(xmin, xmax, &tick1, &dtick, &ndec);

/* X-axis tick marks and labels */

    for (t=tick1; t<xmax; t+=dtick) {
      x1 = conv_coord_x(windata+pos, t);
      XDrawLine(disp, windata[pos].window, windata[pos].gc,  /* Tick mark */
		x1, windata[pos].height-bordy+3,
		x1, windata[pos].height-bordy-3);
      if (fabs(t)< fabs(dtick)/2)
	strcpy(s, "0");
      else{
	if (ndec <= 4 && fabs(t) < 1.e4)
	  sprintf(s, "%-.*f", ndec, t);
	else
	  sprintf(s, "%-.3e", t);
      }

      l=strlen(s); 
      x1 = x1 - XTextWidth(font_info, s, l) / 2 + 1;
      XDrawString(disp, windata[pos].window, windata[pos].gc, x1, 
		  windata[pos].height-bordy+5+font_info->ascent, 
		  s, l);
    }

    X_winborder(winid, 0, 0, 0, 0, FALSE);

/* Restore active color */
    
    if (X_Color_Server) 
      XSetForeground(disp, windata[pos].gc, gcv.foreground);

    XFlush(disp);
  }

  return(OK);
}


/***************************************
 * X title functions 
 **************************************/

static int
X_title(winid_type winid, char *s)
{
  int pos, w, h, x, y;
  XGCValues gcv;
  
  if (!X_is_running)
    return ERROR;
  if (winid_exists(winid, &pos)) {
    if (X_Color_Server) {
      XGetGCValues(disp, windata[pos].gc, GCForeground, &gcv);
      XSetForeground(disp, windata[pos].gc, black);
    }
    XSetClipMask(disp, windata[pos].gc, None);
    w = XTextWidth(font_info, s, strlen(s));
    h = (font_info->ascent + font_info->descent);
    x = (windata[pos].width -windata[pos].xpadl - windata[pos].xpadr - w) / 2
      + windata[pos].xpadl;
    y = h+10;

#ifndef NO_XRECTANGLE
    XClearArea(disp, windata[pos].window, x-10, y, w+20, h*2, FALSE);
    XDrawRectangle(disp, windata[pos].window, 
		   windata[pos].gc, x-10, y, w+20, h*2);
    XDrawRectangle(disp, windata[pos].window, 
		   windata[pos].gc, x-8, y+2, w+16, h*2-4);
    XDrawString(disp, windata[pos].window, windata[pos].gc, x, 
		y+1.5*h-1, s, strlen(s));
#endif

    X_winborder(winid, 0, 0, 0, 0, FALSE);

    /* Restore active color */
    
    if (X_Color_Server) 
      XSetForeground(disp, windata[pos].gc, gcv.foreground);

    XFlush(disp);
    return OK;
  }
  return ERROR;
}



/***************************************
 * X ps functions 
 **************************************/

#define FONTSCALING   1.3

static int
X_ps_setup(winid_type winid, FILE *fp, real *scale, int has_scale, int col_flag) 
{
  real pic_height, xorig, page_width, page_height;
  int pos, xpadl, xpadr, ypadt, ypadb;

  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {
    page_width = 550.0;
    page_height = 775.0;
    xorig = 25.0;

    *scale = page_width/(windata[pos].xmax-windata[pos].xmin);	 
    
    pic_height = (windata[pos].ymax-windata[pos].ymin)* *scale;

/* don't print any decimals in the bounding box */
    fprintf(fp, "%%%%BoundingBox: %.0f %.0f %.0f %.0f\n\n",
	    xorig, (page_height-pic_height)/2, 
	    xorig + page_width, (page_height+pic_height)/2);
    

    if (has_scale == TRUE) {
      get_scale_pad(pos, &xpadl, &xpadr, &ypadt, &ypadb);
      xorig += xpadl*FONTSCALING;
      page_width -= xpadl*FONTSCALING;
      page_height -= ypadb*FONTSCALING;
      
      *scale = page_width/(windata[pos].xmax-windata[pos].xmin);	 
      pic_height = (windata[pos].ymax-windata[pos].ymin)* *scale;
    }
    
    
    /* Define justshow macro (a version of PS 'show' that
       takes the same justifaction parameters, as PL_plopage_string) */
    
    fprintf(fp, "\n%%  /justshow macro: show that takes justification parms\n\n");
    fprintf(fp, "/justshow { dup true charpath flattenpath pathbbox\n");
    fprintf(fp, "3 2 roll sub 3 1 roll sub 5 4 roll -2 div 1.5 add mul\n");
    fprintf(fp, "2 1 roll 4 3 roll 1 add -2 div mul\n");
    fprintf(fp, "currentpoint 3 2 roll add 3 1 roll add 2 1 roll\n");
    fprintf(fp, "moveto show\n");
    fprintf(fp, "} def\n\n");


/* Define marker macros */

    fprintf(fp, "\n%%  Marker macros,  syntax: x y size <cmd>\n\n");
   
    fprintf(fp, "/MakeBox { newpath dup dup dup dup 7 4  roll add\n");
    fprintf(fp, "6 1 roll add 5 4 roll moveto -2 mul 0 rlineto -2 mul\n");
    fprintf(fp, "0 2 1 roll rlineto 2 mul 0 rlineto\n");
    fprintf(fp, "closepath } def\n\n");
    fprintf(fp, "/CBox { MakeBox stroke } def\n\n");
    if (col_flag)
      fprintf(fp, "/CFillBox {  MakeBox fill } def\n\n");
    else
      fprintf(fp, "/CFillBox {gsave 0 setgray MakeBox fill grestore} def\n\n");
    fprintf(fp, "/Asterisk { newpath 3 1 roll moveto dup dup dup dup dup\n");
    fprintf(fp, "0 rmoveto -2 mul 0 rlineto rmoveto -2 mul 0 2 1 roll\n");
    fprintf(fp, "rlineto .7071 mul dup dup dup dup .4142 mul rmoveto\n");
    fprintf(fp, "-2 mul dup -1 mul rlineto 2 mul 0 rmoveto -2 mul\n");
    fprintf(fp, "dup rlineto stroke } def\n\n");
    fprintf(fp, "/Cross {  newpath 3 1 roll moveto \n");
    fprintf(fp, ".7071 mul dup dup dup dup rmoveto -2 mul dup rlineto \n");
    fprintf(fp, "2 mul 0 rmoveto -2 mul dup -1 mul rlineto stroke } def\n\n");
    fprintf(fp, "/Circle { newpath 0 360 arc stroke } def\n\n");
    if (col_flag)
	fprintf(fp, "/FillCircle { newpath 0 360 arc fill } def\n\n");
      else {
	fprintf(fp, "/FillCircle { gsave 0 setgray newpath 0 360 arc fill\n");
	fprintf(fp, "grestore } def\n\n");
      }
    fprintf(fp, "/MakeTriangle { newpath 3 1 roll moveto dup 0 2 1 roll\n");
    fprintf(fp, "rmoveto 0.866 mul dup dup -2 mul rlineto -2 mul 0 rlineto\n");
    fprintf(fp, "closepath } def\n\n");
    fprintf(fp, "/Triangle { MakeTriangle stroke } def\n");
    if (col_flag)
      fprintf(fp, "/FillTriangle { MakeTriangle fill } def\n\n");
    else {
      fprintf(fp, "/FillTriangle { gsave 0 setgray MakeTriangle\n");
      fprintf(fp, "fill grestore } def\n\n");
    }
    fprintf(fp, "/MakeInvTriangle { newpath 3 1 roll moveto dup 0 2 1 roll\n");
    fprintf(fp, "-1 mul rmoveto -0.866 mul dup dup -2 mul rlineto -2 mul\n");
    fprintf(fp, "0 rlineto closepath } def\n\n");
    fprintf(fp, "/InvTriangle { MakeInvTriangle stroke } def\n");
    if (col_flag)
      fprintf(fp, "/FillInvTriangle { MakeInvTriangle fill } def\n\n");
    else {
      fprintf(fp, "/FillInvTriangle { gsave 0 setgray MakeInvTriangle\n");
      fprintf(fp, "fill grestore } def\n\n");
    }
    fprintf(fp, "0.2 setlinewidth\n\n");

/* Load a font */

    fprintf(fp, "/Courier findfont 12 scalefont setfont\n");


/* Translate coordinate system */

    fprintf(fp, "%f %f translate\n", 
	    -windata[pos].xmin* *scale+xorig, 
	    -windata[pos].ymin* *scale+(775-(windata[pos].ymax-windata[pos].ymin)* *scale)/2);
	    
    return OK;
  }

  return ERROR;
}


#define SCALE_CUT 15

static int
X_ps_scale(winid_type winid, FILE *fp, real scale, int col_flag)
{
  real x_orig, y_orig;
  int pos, bordx, dummy, bordy;
  real tick1, dtick, t;
  real xmin, xmax, ymin, ymax;
  int ndec;
  char s[80];

  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {
    fprintf(fp, "\n%% Draw scales\n\n");

    if (col_flag)
      fprintf(fp, "0 0 0 setrgbcolor\n\n");

    x_orig = windata[pos].xmin*scale - SCALE_CUT;
    y_orig = windata[pos].ymin*scale - SCALE_CUT;

    get_scale_pad(pos, &bordx, &dummy, &dummy, &bordy);

/* Coordinate axis */

    fprintf(fp, "newpath\n%f %f moveto\n",
	    x_orig, windata[pos].ymax*scale-SCALE_CUT);
    fprintf(fp, "%f %f lineto\n",
	    x_orig-3, windata[pos].ymax*scale-SCALE_CUT-3);
    fprintf(fp, "%f %f moveto\n",
	    x_orig+3, windata[pos].ymax*scale-SCALE_CUT-3);
    fprintf(fp, "%f %f lineto\n",
	    x_orig, windata[pos].ymax*scale-SCALE_CUT);
    fprintf(fp, "%f %f lineto\n", x_orig, y_orig);
    fprintf(fp, "%f %f lineto\n",
	    windata[pos].xmax*scale-SCALE_CUT, y_orig);
    fprintf(fp, "%f %f lineto\n",
	    windata[pos].xmax*scale-SCALE_CUT-3, y_orig-3);
    fprintf(fp, "%f %f moveto\n",
	    windata[pos].xmax*scale-SCALE_CUT-3, y_orig+3);
    fprintf(fp, "%f %f lineto\n",
	    windata[pos].xmax*scale-SCALE_CUT, y_orig);
    fprintf(fp, "stroke\n");

    fprintf(fp, "newpath\n");
    fprintf(fp, "%f %f moveto\n0 -1 (y) justshow\n",
	    x_orig, windata[pos].ymax*scale-5);

    fprintf(fp, "newpath\n");
    fprintf(fp, "%f %f moveto\n1 0 (x) justshow\n",
	    windata[pos].xmax*scale-5, y_orig);

/* Y-axis tick marks and labels */

    ymin = y_orig / scale;
    ymax = windata[pos].ymax-SCALE_CUT/scale-5/scale;

    tick_dist(ymin, ymax, &tick1, &dtick, &ndec);

    for (t=tick1; t<ymax; t+=dtick) {
      
      fprintf(fp, "newpath\n%f %f moveto\n", x_orig-3, t*scale);
      fprintf(fp, "%f %f lineto\nstroke\n", x_orig+3, t*scale);

      if (fabs(t)< fabs(dtick)/2)
	strcpy(s, "0");
      else{
	if (ndec <= 4 && fabs(t) < 1.e4)
	  sprintf(s, "%-.*f", ndec, t);
	else
	  sprintf(s, "%-.3e", t);
      }

      fprintf(fp, "newpath\n%f %f moveto\n 1 0 (%s) justshow\n",
	      x_orig+10.0-(real) bordx *FONTSCALING, t*scale, s);
    }

    xmin = x_orig / scale;
    xmax = windata[pos].xmax-SCALE_CUT/scale-5/scale;

    tick_dist(xmin, xmax, &tick1, &dtick, &ndec);

/* X-axis tick marks and labels */

    for (t=tick1; t<xmax; t+=dtick) {
      fprintf(fp, "newpath\n%f %f moveto\n", t*scale, y_orig-3);
      fprintf(fp, "%f %f lineto\nstroke\n", t*scale, y_orig+3);

      if (fabs(t)< fabs(dtick)/2)
	strcpy(s, "0");
      else{
	if (ndec <= 4 && fabs(t) < 1.e4)
	  sprintf(s, "%-.*f", ndec, t);
	else
	  sprintf(s, "%-.3e", t);
      }

      fprintf(fp, "newpath\n%f %f moveto\n 0 1 (%s) justshow\n",
	      t*scale, y_orig-10.0, s);
    }
    
    return OK;
  }

  return ERROR;
}

static int
X_ps_data(winid_type winid, FILE *fp, real scale, int col_flag)
{
  int i, pos, ci = -1, line_width = -1, dash = -1;

  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {

    /* Set up clip rectangle */

    fprintf(fp, "\n%% Set up clip rectangle\n\n");

    fprintf(fp, "newpath\n%f %f moveto\n%f %f lineto\n%f %f lineto\n%f %f lineto\nclosepath\nclip\n",
	    windata[pos].xmin* scale, 
	    windata[pos].ymin* scale, 
	    windata[pos].xmin* scale, 
	    windata[pos].ymax* scale, 
	    windata[pos].xmax* scale, 
	    windata[pos].ymax* scale, 
	    windata[pos].xmax* scale, 
	    windata[pos].ymin* scale);

    /* Draw all the lines */

    fprintf(fp, "\n%% Lines\n\n");

    fprintf(fp, "newpath\n");

    for (i=0; i<windata[pos].ls_count; i++) {
      if (ci != windata[pos].ls[i].cind && col_flag) {
	ci = windata[pos].ls[i].cind;
	fprintf(fp, "stroke\n");
	fprintf(fp, "%f %f %f setrgbcolor\nnewpath\n", 
		(real) cm[ci].red / 65535.0,
		(real) cm[ci].green / 65535.0,
		(real) cm[ci].blue / 65535.0);
      }
      if (line_width != windata[pos].ls[i].lw || dash != windata[pos].ls[i].dash) {
	line_width = windata[pos].ls[i].lw;
	dash = windata[pos].ls[i].dash;
	fprintf(fp, "stroke\n");
	fprintf(fp, "%f setlinewidth\n", 0.2 + 0.8*(line_width-1));
	if (dash)
	  fprintf(fp, "[%f] %f setdash\n", 5.0, 0.0); /* turn dashing on */
	else
	  fprintf(fp, "[] %f setdash\n", 0.0); /* turn dashing off */
	fprintf(fp, "newpath\n");
      }
      fprintf(fp, "%f %f moveto\n%f %f lineto\n",
	      windata[pos].ls[i].x0*scale,
	      windata[pos].ls[i].y0*scale,
	      windata[pos].ls[i].x1*scale,
	      windata[pos].ls[i].y1*scale);
    }
    fprintf(fp, "stroke\n");

    
    /* Draw all the text bits */

    fprintf(fp, "\n%% Text\n\n");
    ci = -1;
    
    for (i=0; i<windata[pos].ts_count; i++) {
      if (ci != windata[pos].ts[i].cind && col_flag) {
	ci = windata[pos].ts[i].cind;
	fprintf(fp, "%f %f %f setrgbcolor\n", 
		(real) cm[ci].red / 65535.0,
		(real) cm[ci].green / 65535.0,
		(real) cm[ci].blue / 65535.0);
      }
	fprintf(fp, "newpath\n%f %f moveto\n%d %d (%s) justshow\n",
		windata[pos].ts[i].x*scale,
		windata[pos].ts[i].y*scale,
		windata[pos].ts[i].hjust,
		windata[pos].ts[i].vjust,
		windata[pos].ts[i].s);
    }
    

    /* Draw all the markers */

    fprintf(fp, "\n%% Markers\n\n");

    ci=-1;

    for (i=0; i<windata[pos].ms_count; i++) {
      if (ci != windata[pos].ms[i].cind && col_flag) {
	ci = windata[pos].ms[i].cind;
	fprintf(fp, "%f %f %f setrgbcolor\n", 
		(real) cm[ci].red / 65535.0,
		(real) cm[ci].green / 65535.0,
		(real) cm[ci].blue / 65535.0);
      }
      fprintf(fp, "%f %f ",
	      windata[pos].ms[i].x*scale,
	      windata[pos].ms[i].y*scale);
      switch(windata[pos].ms[i].type) {
      case DOT: 
	fprintf(fp, "0.5 FillCircle\n"); break;
      case CIRCLE: 
	fprintf(fp, "3 Circle\n"); break;
      case BOX: 
	fprintf(fp, "3 CBox\n"); break;
      case FILLED_CIRCLE: 
	fprintf(fp, "3.2 FillCircle\n"); break;
      case FILLED_BOX: 
	fprintf(fp, "3.2 CFillBox\n"); break;
      case ASTERISK: 
	fprintf(fp, "3 Asterisk\n"); break;
      case CROSS: 
	fprintf(fp, "3 Cross\n"); break;
      case TRIANGLE: 
	fprintf(fp, "4 Triangle\n"); break;
      case FILLED_TRIANGLE: 
	fprintf(fp, "4.2 FillTriangle\n"); break;
      case INV_TRIANGLE: 
	fprintf(fp, "4 InvTriangle\n"); break;
      case FILLED_INV_TRIANGLE: 
	fprintf(fp, "4.2 FillInvTriangle\n"); break;
      default: 
	fprintf(stderr, "warning: PS output could not generate marker type %d\n", windata[pos].ms[i].type); break;
      }
    }

    return(OK);
  }
  
  return(ERROR);
}

static int
X_ps_title(winid_type winid, FILE *fp, char *s, real scale, int col_flag)
{
  int pos;

  if (!X_is_running)
    return ERROR;

  if (winid_exists(winid, &pos)) {

    fprintf(fp, "\n%% Title box\n\n");

    fprintf(fp, "%% /Box macro\n\n");
    fprintf(fp, "/Box { dup dup dup 8 7 roll 2 1 roll sub 4 1 roll\n");
    fprintf(fp, "7 6 roll 2 1 roll sub 3 1 roll 6 5 roll add 2 1 roll\n");
    fprintf(fp, " 5 4 roll add 4 copy dup 3 1 roll\n");
    fprintf(fp, "moveto pop pop 2 1 roll lineto\n");
    fprintf(fp, "3 1 roll 3 copy lineto 2 1 roll\n");
    fprintf(fp, "lineto pop pop\n");
    fprintf(fp, " closepath} def\n\n");

    fprintf(fp, "%% /titleshow macro\n\n");
    fprintf(fp, "/titleshow {\n");
    fprintf(fp, "dup true charpath flattenpath pathbbox 4 copy 9 4 roll\n");
    fprintf(fp, "3 2 roll sub 3 1 roll sub dup 4 1 roll 1.5 mul\n");
    fprintf(fp, "2 1 roll dup 4 1 roll -0.5 mul\n");
    fprintf(fp, "currentpoint 3 2 roll add 3 1 roll add 2 1 roll\n");
    fprintf(fp, "9 3 roll 5 4 roll 2 1 roll sub 4 1 roll 2 div dup\n");
    fprintf(fp, "4 3 roll add 10 add 2 1 roll 5 4 roll add 10 sub\n");
    fprintf(fp, "4 1 roll 2 1 roll 4 copy 4 copy 1 setgray\n");
    fprintf(fp, "newpath 2 Box fill stroke 0 setgray\n");
    fprintf(fp, "newpath 2 Box stroke newpath 0 Box stroke\n");
    fprintf(fp, "moveto show } def\n\n");

    if (col_flag)
      fprintf(fp, "0 0 0 setrgbcolor\n\n");

    fprintf(fp, "newpath %f %f moveto (%s) titleshow\n",
	    scale*(windata[pos].xmax+windata[pos].xmin)/2,
	    scale*windata[pos].ymax - 
	    (font_info->ascent + font_info->descent) - 10, s);

    return(OK);
  }
  
  return(ERROR);
}

  
/***************************** E N D  of  X I N T E R **********************/

/*****************************************************************
 *
 * Start of PL routines 
 *
 ****************************************************************/

/*static color_type 
  col_black = {0.0, 0.0, 0.0},
  col_white = {1.0, 1.0, 1.0};
*/

static int PL_is_running = FALSE;
static int n_win;
static winid_type *win;
static real *xm, *ym, *xx, *yx;
static int *scaled;
static char **title;



void PL_start_plot(int w)
{
  if (PL_is_running)
    X_activate(win[w]);
}


void PL_end_plot(void)
{
}



void PL_scale( int w, int yesno)
{
  if (PL_is_running)  {
    if (scaled[w] == TRUE && yesno == FALSE) { /* Switch off scaling */
      X_scale_active(win[w], FALSE);
      X_rescale(win[w], 0, 0, 0, 0, FALSE, same_scale[w]);
      X_redraw(win[w]);
      if (title[w] != NULL)
	X_title(win[w], title[w]);
    }
    if (scaled[w] == FALSE && yesno == TRUE) { /* Switch on scaling */
      X_scale_active(win[w], TRUE);
      X_rescale(win[w], 0, 0, 0, 0, FALSE, same_scale[w]);
      X_redraw(win[w]);
      X_draw_scale(win[w]);
      if (title[w] != NULL)
	X_title(win[w], title[w]);
      
    }
    scaled[w] = (yesno == TRUE) ? TRUE : FALSE;
  }
}


void PL_title( int w, char *s)
{
  if (PL_is_running)  {
    if (s != NULL) {
      X_title(win[w], s);
      if (title[w] == NULL || strcmp(title[w], s)) {
	if (title[w] != NULL) free(title[w]);
	title[w] = (char *) malloc(strlen(s)+1);
	strcpy(title[w], s);
      }
    } else {
      title[w] = NULL;
      X_redraw(win[w]);
      if (scaled[w] == TRUE)
	X_draw_scale(win[w]);
    }
    
  }
}

void PL_window(int w, real xyab[2][2])
{
  if (PL_is_running) {
    X_rescale(win[w], xyab[0][0], xyab[1][0], xyab[0][1], xyab[1][1], TRUE, 
	      same_scale[w]);
    if (scaled[w] == TRUE)
      X_draw_scale(win[w]);
    if (title[w] != NULL)
      X_title(win[w], title[w]);
      
    xm[w] = xyab[0][0];
    ym[w] = xyab[1][0];
    xx[w] = xyab[0][1];
    yx[w] = xyab[1][1];
  }
}

void PL_color(color_type col)
{
  X_color(col);
}

void
PL_start_dash(void){
  windata[active_pos].dash = TRUE;
  return;
}

void
PL_stop_dash(void){
  windata[active_pos].dash = FALSE;
  return;
}

void
PL_line_width(unsigned int line_width){
  if (!X_is_running)
    return;

  XSetLineAttributes(disp, windata[active_pos].gc, line_width, LineSolid, CapButt,
		     JoinMiter);
  windata[active_pos].line_width = line_width;
}

void PL_move(real x_r,real y_r){

  X_moveto(x_r, y_r);
}


void PL_draw(real x_r,real y_r){
  X_lineto(x_r, y_r);
}



void PL_plot_string(char *text, real x, real y, int h_center, int v_center, 
                 int font_number){

/*   write text to plot, note differences between tektronix and postscript
     text                text string to be plotted
     x                   x coordinate where the string should be plotted
     y                   y coordinate where the string should be plotted
     h_center           -1: end the text at x
                         0: center the text around x
                         1: begin the text at x
     v_center           -1: draw the string below y
                         0: draw the string centered around y
                         1: draw the string above y
     font_number         size of font  (currently not used)        */

  X_moveto(x, y);
  X_text(text, h_center, v_center);
}
      
void PL_marker(real x, real y, int type)
{
  X_moveto(x, y);
  X_marker(type);
}


void PL_erase(int w){
/* erase the whole screen */

  if (PL_is_running) {
    X_clear(win[w]); 
    if (title[w] != NULL) free(title[w]);
    title[w] = NULL;
    if (scaled[w] == TRUE)
      X_draw_scale(win[w]);
  }
}

void PL_reset(input_output *io_ptr){
  int i;
  
  /* Shut down X server */

  X_exit(); 

  if (PL_is_running) {
    free(win);
    free(xm);
    free(ym);
    free(xx);
    free(yx);
    free(scaled);
    for (i=0; i<n_win; i++) {
      if (title[i] != NULL) free(title[i]);
    }
    free(title);
  }
  
  PL_is_running = FALSE;

/* close open files */
  if (io_ptr->read_command != stdin){
    fclose(io_ptr->read_command);
    io_ptr->read_command = stdin;
  }

  if (io_ptr->copy_command != NULL){
    fclose(io_ptr->copy_command);
    io_ptr->copy_command = NULL;
/* free the file name */
    free( io_ptr->copy_file_name );
  }

}


void PL_set_same_scale( int window_number, int flag ){
/* check the parameters */
  if (window_number >= 0 && window_number < n_windows && (flag == 0 || flag == 1)){
    same_scale[window_number] = flag;
  }
}


void PL_set_labels( int window_number, char *new_x_label, char *new_y_label ){
/* check the parameters */
  if (window_number >= 0 && window_number < n_windows){
    free( x_label[window_number] );
    free( y_label[window_number] );
    x_label[window_number] = (char *) malloc(sizeof(char)*(strlen(new_x_label)+1));
    y_label[window_number] = (char *) malloc(sizeof(char)*(strlen(new_y_label)+1));
    strcpy( x_label[window_number], new_x_label );
    strcpy( y_label[window_number], new_y_label );
  }
}


int PL_init(input_output *io_ptr, char **titles, real xy[][2], real *w, 
	    real *h, int n, int use_graphics){
  int i;


  n_win =n;

  /* Initially read the commands from standard in and don't save the commands */
  
  io_ptr->read_command = stdin;
  io_ptr->copy_command = NULL;
  
  /* Make sure to shut down any pre existing windows */
  
  X_exit(); 
  if (PL_is_running) {
    free(win);
    free(xm);
    free(ym);
    free(xx);
    free(yx);
    free(scaled);
    free(same_scale);
    for (i=0; i<n; i++) {
      if (title[i] != NULL) free(title[i]);
    }
    free(title);
    free(x_label);
    free(y_label);
  }

  PL_is_running = FALSE;

  /* Check if the application doesn't need any windows */

  if (n_win == 0 || !use_graphics) 
    return ERROR;

  /* Start up X */
  
  if (X_init() == ERROR){
    printf("Warning: No graphical output because the X-server is not responding\n");
    return ERROR;
  }

  PL_is_running = TRUE;

/* Allocate window data memory */

  win     = (winid_type *) malloc(sizeof(winid_type)*n);
  xm      = (real *) malloc(sizeof(real)*n);
  ym      = (real *) malloc(sizeof(real)*n);
  xx      = (real *) malloc(sizeof(real)*n);
  yx      = (real *) malloc(sizeof(real)*n);
  scaled  = (int *) malloc(sizeof(int)*n);
  title   = (char **) malloc(sizeof(long)*n);
  same_scale = (int *) malloc(sizeof(int)*n);
  x_label = (char **) malloc(sizeof(char *)*n);
  y_label = (char **) malloc(sizeof(char *)*n);


/* Open windows */

  for (i=0; i<n_win; i++) {
    same_scale[i] = 1;
    x_label[i] = (char *) malloc(2*sizeof(char));
    y_label[i] = (char *) malloc(2*sizeof(char));
    strcpy( x_label[i], "x" );
    strcpy( y_label[i], "y" );
    win[i] = X_winopen(titles[i], xy[i][0], xy[i][1], w[i], h[i]);
    xm[i] = -1.0;
    ym[i] = -1.0;
    xx[i] =  1.0;
    yx[i] =  1.0;
    scaled[i] = FALSE;
    title[i] = NULL;
  }
  
  return OK;
}


void 
PL_poll(void)
{
  winid_type w;
  real x0, y0, x1, y1, object_size;
  int ret, wint=0, i;

  if (!PL_is_running) 
    return;

  if (X_poll_resize(&w) == TRUE) {
    for (i=0; i<n_win; i++)
      if (win[i]==w) { 
	wint = i;
	i = n_win;
      }

    X_rescale(w, 0.0, 0.0, 0.0, 0.0, FALSE, same_scale[wint]);
    X_winborder(w, 0, 0, 0, 0, FALSE);
    X_redraw(w);
    if (scaled[wint] == TRUE)
      X_draw_scale(w);
    if (title[wint] != NULL)
      X_title(w, title[wint]);
  }
  if ((ret=X_poll_zoom(&w, &x0, &y0, &x1, &y1)) != FALSE) {
    for (i=0; i<n_win; i++) 
      if (win[i]==w) { 
	wint = i;
	i = n_win;
      }

/* try to make sure that the user can't zoom in or out indefinitely */
    object_size = real_max( xx[wint] - xm[wint], yx[wint] - ym[wint] );
    if (ret == TRUE && 
	x1-x0 > 1.e-5*object_size &&
	y1-y0 > 1.e-5*object_size &&
	x1-x0 < 10.0*object_size &&
	y1-y0 < 10.0*object_size ) {
      X_rescale(w, x0, y0, x1, y1, TRUE, same_scale[wint]);
      X_winborder(w, 0, 0, 0, 0, FALSE);
      X_redraw(w);
      if (scaled[wint] == TRUE)
	X_draw_scale(w);
      if (title[wint] != NULL)
	X_title(w, title[wint]);
    }
    if (ret == RESET) {
      X_rescale(w, xm[wint], ym[wint], xx[wint], yx[wint], TRUE, same_scale[wint]);
      X_winborder(w, 0, 0, 0, 0, FALSE);
      X_redraw(w);
      if (scaled[wint] == TRUE)
	X_draw_scale(w);
      if (title[wint] != NULL)
	X_title(w, title[wint]);
    }
  }
}


void
PL_postscript(int w, FILE *fp, int col_flag)
{
  real sc;

  if (PL_is_running) {
    fprintf(fp, "%%!\n");
    fprintf(fp, "%%%%Creator: PL_postscript\n");
    fprintf(fp, "%%%%Title: PL window dump\n\n");
    X_ps_setup(win[w], fp, &sc, scaled[w], col_flag);
    if (scaled[w] == TRUE)
      X_ps_scale(win[w], fp, sc, col_flag);

    X_ps_data(win[w], fp, sc, col_flag);
    if (title[w] != NULL) 
      X_ps_title(win[w], fp, title[w], sc, col_flag);

    fprintf(fp, "showpage\n");
  }
}



/***************************** E N D  of  P L **********************/

/***************************************************************************
 *
 * Start of command interpreter 
 *
 **************************************************************************/



/* private functions. make them static so they can't be used outside this file */
static char *gettok(input_output *io_ptr, char *prompt,
                    char **cmd, int n);
static char *getline_no_com( char *s, int n, input_output *io_ptr );


int 
get_command(input_output *io_ptr, char *prompt,
	    char **command, int n_command,
	    int command_level, int *save_on_copy, int *argument) {
  char *token, *cnull, *first_token, *second_token, *second_command=NULL;
  char *third_token=NULL, *third_command=NULL;
  int match, unique, i, first_length, second_length, third_length;
  input_output pause_io;
  int pause_icom, quit;
  char *pause_com[3], *pause_help[3];

  match = -1;
  unique = 1;

  token = gettok( io_ptr, prompt, command, n_command + 1);

  /* check for "?" */
  
  while (token[0] == '?' || token[0] == '\n' ) {
    if (token[0] == '?'){
      printf("Choose one of the following commands:\n");
      for (i=0; i<=n_command; i++){
	printf("%s\n",command[i]);
	if ((i+1)%3 == 0) printf("\n");
      }
    }
/* read a new answer */
    token = gettok( io_ptr, prompt, command, n_command + 1 );
  }

/* check for empty command */

  if (token[0] == '\n') 
    return -1;
  else{      /* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* check for a pause */
  if (strcmp(token, "pause") == 0){
    pause_io.read_command = stdin;
    pause_io.copy_command = NULL;
    pause_com[0] = "proceed";
    pause_com[1] = "break";
    pause_com[2] = "help";
    pause_help[0] = "Continue with the tutorial, i.e. read next instruction from "
      "the command file.";
    pause_help[1] = "Terminate the tutorial session, i.e. close the command file "
      "and start reading commands from standard input.";
    pause_help[2] = NULL;
    quit = 0;
    do{
      pause_icom = get_command( &pause_io, "Pause>", pause_com, 2, 
			       0, NULL, NULL);

      switch (pause_icom){

      case 0:
/* proceed */
	quit = 1;
	break;

      case 1:
/* break */
	fclose(io_ptr->read_command);
	io_ptr->read_command = stdin;
	quit = 1;
	break;

      case 2:
/* help */
	while ((pause_icom = get_command( &pause_io, 
					 "help on subject (? lists all subjects)>", 
					 pause_com, 2, 0, NULL, NULL)) == -1);
	if (pause_help[pause_icom] == NULL)
	  general_help(); 
	else
	  print_help( pause_help[pause_icom] ); 
	break;

      default:
	;
      }

    }
    while( !quit );
    return -1;
  }	 
  
  /* decompose the token into its sub parts */
    
  first_token = token;
  first_length = strcspn( token, "-" );

  if ( (second_token = strchr( token, '-' )) == NULL ){
    second_length = 0;
    third_length = 0;
  } else {
    second_length = strcspn( &(second_token[1]), "-" ) + 1;
    if ( (third_token = strchr( &(second_token[1]), '-' )) == NULL ){
      third_length = 0;
    }
    else{
      third_length = strcspn( &(third_token[1]), "-" ) + 1;
    }
  }
  
/* check for matches */

  for (i=0; i<=n_command; i++){
/* new algorithm */
/* decompose the command */
    if (second_length > 0){
      second_command = strchr( command[i], '-' );
      third_command = NULL;
      if ( third_length > 0){
	if (second_command == NULL)
	  third_command = NULL;
	else
	  third_command = strchr( &(second_command[1]), '-' );
      }
    }
/* the first part of the command matches, then check if either the second part
 is absent or if it also matches */
    if ( strncmp( first_token, command[i], first_length)==0 &&
	(second_length == 0 ||
	 (second_command != NULL && strncmp( second_token, second_command, 
					    second_length )==0)) &&
	(third_length == 0 ||
	 (third_command != NULL && strncmp( third_token, third_command, 
					   third_length )==0)) ){
	if (match == -1)  
	  unique = 1; 
	else 
	  unique = 0; 
	match = i; 
      }
/* old algorithm */
/*     length = strlen( token ); */
/*     if ( strncmp( token, command[i], length)==0 ){ */
/*       if (match == -1)  */
/* 	unique = 1; */
/*       else */
/* 	unique = 0; */
/*       match = i; */
  }
/* even if the command wasn't unique, it may be a perfect match. This happens
   for instance when commmand[0]="s", command[1]="show" and token="s". */
  if (match != -1 && !unique){
    match = -1;
    for (i=0; i<=n_command; i++)
      if (strcmp( token, command[i] )==0){
	match=i;
      }
  }


  if (match==-1){
    if (!unique)
      printf("Not an unique command: %s\n", token);
    else
      printf("Unknown command: %s\n", token);
  }

/* copy the command onto the copy stream */
  if (match != -1 && io_ptr->copy_command != NULL &&
      (save_on_copy==NULL || save_on_copy[match])){

/* indent accordning to command level */
    if (strcmp("exit", command[match])==0)
      command_level += -1;
    for (i=0; i< command_level; i++)
      fprintf( io_ptr->copy_command, "\t");

/* output the command */
    fprintf( io_ptr->copy_command, "%s", command[match] );
/* can we expect an argument to this command which should be printed on the
   same line? */
    if (argument == NULL || !argument[match])
      fprintf( io_ptr->copy_command, "\n");
    else
      fprintf( io_ptr->copy_command, " ");
  }

/* echo a valid command on stdout */
  if (match != -1 ){
/* output the command */
    printf( "%s", command[match] );
    if ( io_ptr->read_command == stdin )
      printf( "\n" );
    else{
/* if the command was read from a file:
   can we expect an argument to this command which should be printed on the
   same line? */
      if (argument == NULL || !argument[match])
	printf( "\n" );
      else
	printf( " " );
    }
  }

  return match;

}


char 
*get_word(input_output *io_ptr, char *prompt, char *deflt,
	  int save_on_copy){
  char *token, *cnull, prompt2[120];

  sprintf(prompt2, "%s (%s)", prompt, deflt );
  do{
    token = gettok( io_ptr, prompt2, NULL, 0 );
  }    
  while ( token[0] == '?' );

  if (token[0] == '\n') 
    token = deflt;
  else{
/* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* copy it onto the copy stream */
  if ( io_ptr->copy_command != NULL && save_on_copy )
    fprintf( io_ptr->copy_command, "%s\n", token );

/* echo the command on stdout if the command was read from file */
  if ( io_ptr->read_command != stdin )
    printf( "%s\n", token );

  return token;
}

int get_int(input_output *io_ptr, char *prompt, int deflt, int level){
  char *token, *cnull, prompt2[120];
  int i;

  sprintf(prompt2, "%s (%i)", prompt, deflt );
  do{
    token = gettok( io_ptr, prompt2, NULL, 0 );
/* check that token consists only of (+/-) and digits. Otherwise return NULL */
  }    
  while ( token[0] == '?' );

  if (token[0] == '\n') 
    sprintf(token,"%i", deflt);
  else{
/* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* copy it onto the copy stream */
  if ( io_ptr->copy_command != NULL ){
/* indent according to the command level */
    for (i=0; i< level; i++)
      fprintf( io_ptr->copy_command, "\t");
    fprintf( io_ptr->copy_command, "%s\n", token );
  }

/* echo the command on stdout if the command was read from file */
  if ( io_ptr->read_command != stdin )
    printf( "%s\n", token );

  return atoi(token);
}

real get_real(input_output *io_ptr, char *prompt, real deflt, int level){
  char *token, *cnull, prompt2[120];
  int i;

  sprintf(prompt2, "%s (%e)", prompt, deflt );
  do{
    token = gettok( io_ptr, prompt2, NULL, 0 );
/* check that token can be converted to a real number. Otherwise return NULL */
  }    
  while ( token[0] == '?' );

  if (token[0] == '\n') 
    sprintf(token,"%e", deflt);
  else{
/* replace any '\n' by '\0' */
    if ( (cnull = strchr(token,'\n')) != NULL )
      cnull[0] = '\0';
  }

/* copy it onto the copy stream */
  if ( io_ptr->copy_command != NULL ){
/* indent according to the command level */
    for (i=0; i< level; i++)
      fprintf( io_ptr->copy_command, "\t");
/* output the number */
    fprintf( io_ptr->copy_command, "%s\n", token );
  }

/* echo the command on stdout if the command was read from file */
  if ( io_ptr->read_command != stdin )
    printf( "%s\n", token );

  return atof(token);
}

static char *gettok(input_output *io_ptr, char *prompt,
		    char **command, int n_command) {
  static char line[MAX_LINE_SIZE];
  char *token;
  static int first_call=TRUE;

/* when strtok is called with line as first parameter, it starts parsing from the
   first character of line and reads until it finds one of the non-tokens that are
   given as the second parameter. When strtok is called with NULL as the first 
   parameter, it continues to read from the previous line (which it stored internally)
   until it finds one of the non-tokens. We want the following functionality from
   this routine:
   when a line only contains a newline, it should be returned as a token.
   when a line contains other characters, newline should not be regarded as a token.
   This is because if the user by misstake gives the line `command <space> <return>',
   we want it to be interpreted only as `command'.
   To accomplish this, we give ' \t\n' as non-tokens when NULL is the first argument
   to strtok, while we only give ' \t' as non-tokens when 'line' is the first 
   argument. */

/* give the graphics a chance to update the plot for size changes */
  PL_poll();

/* strtok shold not be called with a NULL argument the very first time! */
  if (first_call || (token = strtok(NULL," \t\n")) == NULL){
    first_call=FALSE;
    do{
/* output the prompt only if we read from stdin */
      if (io_ptr->read_command == stdin){
	getline_no_block( prompt, line, command, n_command );
      }
      else{
	getline_no_com( line, 120, io_ptr );
	PL_poll();
      }
    }
    while( (token = strtok(line," \t"))==NULL );
  }
  return token;
}

static char *getline_no_com(char *line_no_comment, int buf_size, 
			    input_output *io_ptr ){
  char *first_comment;
  int i;

  if (fgets( line_no_comment, buf_size, io_ptr->read_command ) == NULL){
/* switch to stdin at end-of-file. */
    if (feof(io_ptr->read_command)){
/* If end-of-file occurs on standard in, try to reset it. */
      if (io_ptr->read_command == stdin){
	clearerr( stdin );
	printf("\n");
      }
/* if end-of-file occurs for a command file, close the command file and switch to stdin */
      else{
	fclose(io_ptr->read_command);
	io_ptr->read_command = stdin;
      }
/* important to reset the line */
      line_no_comment[0]='\0';
    }
  }

/* is there a comment ? */
  if ((first_comment=strchr( line_no_comment, '#' ))!=NULL){
/* echo the comment on stdout */
    printf("%s", first_comment);
/* overwrite the first occurance of the comment by '\0' */
/* and fill the rest of the string with blanks */
    for (i=0; first_comment[i] != '\0'; i++)
      first_comment[i] = ' ';
    *first_comment = '\0';
  }

/* return the pointer to the possibly modified string */
  return line_no_comment;
}

 
/********************************** START *********************************/


static void
set_raw(int on){
  struct termio  arg ;
  static char    min ;
  static char    this_time;
  
  if (on != RAWOFF)
    {
      ioctl (fileno(stdin), TCGETA, &arg);
      arg.c_lflag &= ~ICANON & ~ECHO;
      min  = arg.c_cc[VMIN];
      this_time = arg.c_cc[VTIME];
      arg.c_cc[VMIN]  = (on == RAWON) ? 1 : 0;
/* block the read for 0.1 seconds instead of 0.0 seconds which made the program
   a real cpu eater when it was idle */
      arg.c_cc[VTIME] = 1;
      ioctl (fileno(stdin), TCSETA, &arg);
    }
  else
    {
      ioctl (fileno(stdin), TCGETA, &arg);
      arg.c_lflag |= ICANON | ECHO;
      arg.c_cc[VMIN]  = min ;
      arg.c_cc[VTIME] = this_time;
      ioctl (fileno(stdin), TCSETA, &arg);
    }
}

/* Returns how many (initial) charachters of strings s1 and s2 match */
/* i=strlen(s1)=strlen(s2) if s1==s2 */

static int
strolof(char *s1, char *s2) {
  int i=0;
  while (*s1 != 0 && *s2 != 0 && *s1++ == *s2++)
    i++;
  return i;
}


void
getline_no_block( char *prompt, char *buf, char **cmds, int n) {
  int c, i, k, cancel_flag=FALSE, last_esc, cb_pos, new_pos;
  unsigned char ch, bell=7;
  static char backsp[3] = {8, 32, 8};
  char compl[MAX_LINE_SIZE];
  static char cmd_buf[SCROLLBACK_SIZE][MAX_LINE_SIZE];
  static int cmd_buf_pos = 0, cmd_buf_full = FALSE;

  do {

  last_esc = -1;
  printf("%s", prompt);
  fflush(stdout);

  for (i=0; i<n; i++)
    if (!strcmp(cmds[i], "cancel")) {
      cancel_flag = TRUE;
      i=n;
    }

  set_raw(RAWNDELON);
    
  c=0;
  cb_pos = cmd_buf_pos;
  cmd_buf[cb_pos][0] = 0; /* Clear current slot in command buffer */

  do {
    while (read(fileno(stdin), &ch, 1) == 0)
      PL_poll();
    
    /* Trap ESC-sequences generated by arrow keys */

    if (last_esc == 0 && ch == '[') {
      last_esc++;
      ch=0;
    } else {
      if (last_esc == 1 && ch == 'A') ch = 'P' & 31; /* UP   --> CTRL-P */
      if (last_esc == 1 && ch == 'B') ch = 'N' & 31; /* DOWN --> CTRL-N */
      last_esc = -1;
    }
    if (ch == 27)
      last_esc = 0;

    /* Print printable characters */

    if (ch >= 32 && ch < 256 && ch != 127) { /* Should be able to use 8 bit input */
      buf[c++]=ch;
      write(fileno(stdout), &ch, 1);
    }

    /* Backspace or Delete */

    if ((ch == 8 || ch == 127) && (c > 0)) {  
      write(fileno(stdout), backsp, 3);
      c--;
    }

    /* Completion! */

    if (ch == '\t') { 
      k=c-1;
      buf[c]=0;
      strcpy(compl, buf);
      for (i=0; i<n; i++) {
	if (strlen(cmds[i])> c && !strncmp(buf, cmds[i], c)) {
	  if (k==c-1) { /* First matching command */
	    strcpy(compl, cmds[i]);
	    k=strlen(compl);
	  } else { /* Subsequent matching commands */
	    k = strolof(compl, cmds[i]);
	    compl[k]=0;
	  }
	}
      }
      if (k <= c)
	write(fileno(stdout), &bell, 1);
      else {
	write(fileno(stdout), compl+c, k-c);
	strncpy(buf, compl, k);
	c=k;
      }
    }

    /* CTRL-P, previous command */
    
    if (ch == ('P' & 31)) { 
      new_pos = (cb_pos + SCROLLBACK_SIZE - 1) % SCROLLBACK_SIZE;
      if ((new_pos != SCROLLBACK_SIZE-1 || cmd_buf_full == TRUE) &&
	  new_pos != cmd_buf_pos) {
	cb_pos = new_pos;

	/* Erase current line */

	memset(buf, 8, c); write(fileno(stdout), buf, c);
	memset(buf, 32, c); write(fileno(stdout), buf, c);
	memset(buf, 8, c); write(fileno(stdout), buf, c);

	strcpy(buf, cmd_buf[cb_pos]);
	write(fileno(stdout), cmd_buf[cb_pos], strlen(cmd_buf[cb_pos]));
	c = strlen(cmd_buf[cb_pos]);
      } else
	write(fileno(stdout), &bell, 1);   /* Hit end of cmd_buf */
    }

    /* CTRL-N, next command */

    if (ch == ('N' & 31)) { 
      if (cb_pos != cmd_buf_pos) {
	cb_pos = (cb_pos + 1) % SCROLLBACK_SIZE;

	/* Erase current line */

	memset(buf, 8, c); write(fileno(stdout), buf, c);
	memset(buf, 32, c); write(fileno(stdout), buf, c);
	memset(buf, 8, c); write(fileno(stdout), buf, c);

	strcpy(buf, cmd_buf[cb_pos]);
	write(fileno(stdout), cmd_buf[cb_pos], strlen(cmd_buf[cb_pos]));
	c = strlen(cmd_buf[cb_pos]);
      } else
	write(fileno(stdout), &bell, 1);   /* Hit start of cmd_buf */
    }
      
    /* CTRL-D == cancel (if that is among commands in command array) */

    if (ch == ('D' & 31) && cancel_flag == TRUE) {    
      strcpy(buf, "cancel");
      c=6;
      ch='\n';
    }
    
  } while (ch != '\n');
  
  write(fileno(stdout), &ch, 1);  /* Print the last newline */
  
  buf[c++] = '\n';
  buf[c]   = '\0';

  /* Save command in command buffer */

  if (c != 1) { /* Not empty command -> save in command buffer */
    strncpy(cmd_buf[cmd_buf_pos], buf, c-1);
    cmd_buf[cmd_buf_pos][c-1] = 0;
    cmd_buf_pos = (cmd_buf_pos + 1) % SCROLLBACK_SIZE;
    if (cmd_buf_pos == 0)
      cmd_buf_full = TRUE;
  }
    
  set_raw(RAWOFF);

  /* Special! If first character is '!': escape to shell and
     return nothing to (hide it from) command interpreter */

  if (buf[0] == '!')
    system(buf+1);
} while (buf[0] == '!');
}


/********************************** END *********************************/


void
wait_for_key(char *prompt) 
{
  unsigned char ch;

  printf("%s", prompt);
  fflush(stdout);

  set_raw(RAWNDELON);
    
  do {
    while (read(fileno(stdin), &ch, 1) == 0)
      PL_poll();
  } while (ch != '\n' && ch != ' ');

  ch = '\n';

  write(fileno(stdout), &ch, 1);  /* Print newline */

  set_raw(RAWOFF);

}

void get_new_color( color_type *color_ptr, char type ){
  const int n_color=7;
  static int im=0, ic=0;
  int i;

  if (type == 'c'){
    i  = ic;
    ic = (ic+1) % n_color;
  }
  else{
    i  = im;
    im = (im+1) % n_color;
  }

  
  switch( i ){

  case 0:
/* red */
    color_ptr->r = 255.0;
    color_ptr->g = 0.0;
    color_ptr->b = 0.0;
    break;

  case 1:
/* green */
    color_ptr->r = 50.0;
    color_ptr->g = 205.0;
    color_ptr->b = 50.0;
    break;

  case 2:
/* blue */
    color_ptr->r = 0.0;
    color_ptr->g = 0.0;
    color_ptr->b = 255.0;
    break;

  case 3:
/* brown */
    color_ptr->r = 190.0;
    color_ptr->g = 81.0;
    color_ptr->b = 6.0;
    break;

  case 4:
/* purple */
    color_ptr->r = 240.0;
    color_ptr->g = 26.0;
    color_ptr->b = 207.0;
    break;

  case 5:
/* black */
    color_ptr->r = 0.0;
    color_ptr->g = 0.0;
    color_ptr->b = 0.0;
    break;

  case 6:
/* grey */
    color_ptr->r = 155.0;
    color_ptr->g = 166.0;
    color_ptr->b = 187.0;
    break;

  default:
    ;
  }

  color_ptr->r = color_ptr->r/255.0;
  color_ptr->g = color_ptr->g/255.0;
  color_ptr->b = color_ptr->b/255.0;

  return;
}

int get_yes_no(input_output *io_ptr, char *prompt, int level ){
#define ncom 1
  char *command[ncom+1];
  int icom, quit, yes_no;
  int *argument=NULL, *save_on_copy=NULL;

  command[0] ="yes";
  command[1] ="no";

  quit = 0;
  yes_no = 0;
  do{
    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    switch (icom) {
    case 0: 
      yes_no = 1;
      quit = 1;
      break;
    case 1:
      yes_no = 0;
      quit = 1;
      break;
    default:
      ;
    }
  }
  while(!quit);

  return yes_no;
}

int get_on_off(input_output *io_ptr, char *prompt){
#define ncom 1
  char *command[ncom+1];
  int icom, quit, on_off;
  int *argument=NULL, *save_on_copy=NULL, level=-1;

  command[0] ="on";
  command[1] ="off";

  quit = 0;
  on_off = 0;
  do{
    icom = get_command( io_ptr, prompt, command, ncom, 
		       level, save_on_copy, argument);

    switch (icom) {
    case 0: 
      on_off = 1;
      quit = 1;
      break;
    case 1:
      on_off = 0;
      quit = 1;
      break;
    default:
      ;
    }
  }
  while(!quit);

  return on_off;
}
