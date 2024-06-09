/* $Id: db.h,v 1.2 1996/05/23 18:32:30 andersp Exp $ */
/* db.h - header file for db				*/

#include <macros.h>
#include <string.h>
#include "real.h"
#include "c_array.h"

#define NAMLEN 80
#define DBMAX  5
#define MAXRANK  7
#define DEFAULT_CLASS "group"

/* Data */

  typedef struct {
    char name[NAMLEN];
    int  sd_id;
    int  h_id;
    int  vg_ref;
    int  root;
  } group_descriptor;

  group_descriptor *group_list[DBMAX];

/* Function declarations */

int dbopen(int gid, char *path); 
int dbclose(int gid);
int dbput_array(int gid,   char *name, int numtype,
                int rank,  int *dims,  int it,
                void *a);
int dbget_array(int gid,   char *name, int numtype,
                int rank,  int *dims,  int it,
                void *a);
group_descriptor *dbget_group(int gid);

/* c_array interface */
int dbput_int(int gid,char *name,int *a,int it,int itflg);
int dbput_int_array_1d(int gid,char *name,int_array_1d *a,int it,int itflg);
int dbput_int_array_2d(int gid,char *name,int_array_2d *a,int it,int itflg);
int dbput_int_array_3d(int gid,char *name,int_array_3d *a,int it,int itflg);
int dbput_int_array_4d(int gid,char *name,int_array_4d *a,int it,int itflg);

int dbput_real(int gid,char *name,real *a,int it,int itflg);
int dbput_real_array_1d(int gid,char *name,real_array_1d *a,int it,int itflg);
int dbput_real_array_2d(int gid,char *name,real_array_2d *a,int it,int itflg);
int dbput_real_array_3d(int gid,char *name,real_array_3d *a,int it,int itflg);
int dbput_real_array_4d(int gid,char *name,real_array_4d *a,int it,int itflg);

int dbget_int(int gid,char *name,int *a,int it,int itflg);
int dbget_int_array_1d(int gid,char *name,int_array_1d *a,int it,int itflg);
int dbget_int_array_2d(int gid,char *name,int_array_2d *a,int it,int itflg);
int dbget_int_array_3d(int gid,char *name,int_array_3d *a,int it,int itflg);
int dbget_int_array_4d(int gid,char *name,int_array_4d *a,int it,int itflg);

int dbget_real(int gid,char *name,real *a,int it,int itflg);
int dbget_real_array_1d(int gid,char *name,real_array_1d *a,int it,int itflg);
int dbget_real_array_2d(int gid,char *name,real_array_2d *a,int it,int itflg);
int dbget_real_array_3d(int gid,char *name,real_array_3d *a,int it,int itflg);
int dbget_real_array_4d(int gid,char *name,real_array_4d *a,int it,int itflg);
