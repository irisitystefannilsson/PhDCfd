# 1 "hdf_stuff.c"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "hdf_stuff.c"




static int file_id=-1, sd_id=-1;
static char *access_mode;

int32
open_hdf_file(char * name, char access){
  int32 root, ref;
  int fd_tmp;
  char vname[VGNAMELENMAX];

  if (file_id > 0){
    printf("open_hdf_file: ERROR: cannot open an open file!\n");
    return -1;
  }

  if (access == 'i'){

    access_mode = "w";

    if ((sd_id = SDstart(name, DFACC_CREATE)) <= 0){
      return -1;
    }

    file_id = Hopen(name, DFACC_RDWR, 0);
    Vstart(file_id);

    root = Vattach(file_id, -1, access_mode);

    Vsetname(root, "root");
    Vsetclass(root, "directory");
  }
  else if (access == 'r'){

    if ((fd_tmp=open(name, O_RDONLY, 0)) == -1){
      printf("Error: could not open the file `%s'\n", name);
      return -1;
    }
    else{
      close(fd_tmp);
    }
    printf("Opening a HDF data-base...");

    access_mode = "r";

    if ((sd_id = SDstart(name, DFACC_READ)) <= 0){
      return -1;
    }

    file_id = Hopen(name, DFACC_READ, 0);
    Vstart(file_id);


    ref = Vgetid(file_id, -1);
    root = Vattach(file_id, ref, access_mode);
    Vgetname(root,vname);
    if( root <= 0 || strcmp(vname,"root") ){
      printf("\nopen_hdf_file: ERROR: There is no `root' directory in this "
      "dataBase file!\n");

      Vend(file_id);

      Hclose(file_id);

      SDend(sd_id);

      file_id = sd_id = -1;
      return -1;
    }
    printf("done\n");
  }
  else{
    printf("open_hdf_file: ERROR: Unknown access mode %c \n", access);
    return -1;
  }


  return root;
}

int32
get_file_id(void){
  return file_id;
}

int32
get_sd_id(void){
  return sd_id;
}

void
close_hdf_file(int32 root){
  if (file_id == -1){
    printf("close_hdf_file: ERROR: Cannot close a closed file!\n");
    return;
  }

  Vdetach(root);

  Vend(file_id);

  Hclose(file_id);

  SDend(sd_id);


  file_id = sd_id = -1;
}

int32
create_dir(char * name, char * class_name, int vgroup_id){
# 123 "hdf_stuff.c"
  int32 sub_directory;

  if (file_id == -1){
    printf("ERROR: There is no open database file!\n");
    return -1;
  }


  if ( (sub_directory = Vattach(file_id, -1, "w"))<=0 ){
    printf("create: FATAL ERROR in creating a new directory!\n");
    return -1;
  }

  Vsetname(sub_directory, name);
  Vsetclass(sub_directory, class_name);


  Vinsert(vgroup_id, sub_directory);

  return sub_directory;
}

int32
locate_dir(char * name, int vgroup_id){







  int32 tag, ref, sub_directory=-1;
  char vname[VGNAMELENMAX];
  int found=FALSE, npairs, i, status;

  if (file_id == -1){
    printf("ERROR: There is no open database file!\n");
    return -1;
  }


  npairs = Vntagrefs(vgroup_id);
  for(i=0; i<npairs; i++ ){

    status = Vgettagref(vgroup_id, i, &tag, &ref );
    if( Visvg(vgroup_id, ref) ){

      sub_directory = Vattach(file_id, ref, access_mode);
      Vgetname(sub_directory,vname);
      if( !strcmp(name, vname) ){
 found=TRUE;
 break;
      }
      else
        Vdetach(sub_directory);
    }
  }
  if( !found ){
    printf("locate: ERROR: unable to find directory %s\n", name);
  }

  return sub_directory;
}
# 218 "hdf_stuff.c"
int hput_int ( int x, char *name, int vgroup_id ){ int num=1; int32 vdata_id; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } vdata_id = VSattach(file_id, -1, access_mode); VSsetname(vdata_id, name); VSfdefine(vdata_id, name, DFNT_INT32, 1); VSsetfields(vdata_id, name); VSwrite(vdata_id, (unsigned char*)(&x), num, FULL_INTERLACE); Vinsert(vgroup_id, vdata_id); VSdetach(vdata_id); return TRUE; }

int hput_real ( real x, char *name, int vgroup_id ){ int num=1; int32 vdata_id; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } vdata_id = VSattach(file_id, -1, access_mode); VSsetname(vdata_id, name); VSfdefine(vdata_id, name, DFNT_REAL, 1); VSsetfields(vdata_id, name); VSwrite(vdata_id, (unsigned char*)(&x), num, FULL_INTERLACE); Vinsert(vgroup_id, vdata_id); VSdetach(vdata_id); return TRUE; }

int hput_float ( float x, char *name, int vgroup_id ){ int num=1; int32 vdata_id; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } vdata_id = VSattach(file_id, -1, access_mode); VSsetname(vdata_id, name); VSfdefine(vdata_id, name, DFNT_FLOAT32, 1); VSsetfields(vdata_id, name); VSwrite(vdata_id, (unsigned char*)(&x), num, FULL_INTERLACE); Vinsert(vgroup_id, vdata_id); VSdetach(vdata_id); return TRUE; }
int hput_double ( double x, char *name, int vgroup_id ){ int num=1; int32 vdata_id; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } vdata_id = VSattach(file_id, -1, access_mode); VSsetname(vdata_id, name); VSfdefine(vdata_id, name, DFNT_FLOAT64, 1); VSsetfields(vdata_id, name); VSwrite(vdata_id, (unsigned char*)(&x), num, FULL_INTERLACE); Vinsert(vgroup_id, vdata_id); VSdetach(vdata_id); return TRUE; }
# 270 "hdf_stuff.c"
int hget_int ( int *x, char *name, int vgroup_id ) { int npairs, found=FALSE, i, status; int32 vdata_tag, vdata_ref, vdata_id, n_records, nt; char vdata_name[VSNAMELENMAX]; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } npairs = Vntagrefs(vgroup_id); for(i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &vdata_tag, &vdata_ref ); if( Visvs(vgroup_id, vdata_ref) ){ vdata_id=VSattach(file_id, vdata_ref, access_mode); VSgetname(vdata_id, vdata_name); nt = VFfieldtype(vdata_id, 0); if(!strcmp(name, vdata_name) && nt == DFNT_INT32){ found=TRUE; VSQuerycount(vdata_id, &n_records); VSread(vdata_id, (unsigned char*)x, n_records, FULL_INTERLACE ); VSdetach(vdata_id); break; } VSdetach(vdata_id); } } if( !found ){ printf("hget_" "int" ": ERROR searching for %s\n", name); } return found; }

int hget_real ( real *x, char *name, int vgroup_id ) { int npairs, found=FALSE, i, status; int32 vdata_tag, vdata_ref, vdata_id, n_records, nt; char vdata_name[VSNAMELENMAX]; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } npairs = Vntagrefs(vgroup_id); for(i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &vdata_tag, &vdata_ref ); if( Visvs(vgroup_id, vdata_ref) ){ vdata_id=VSattach(file_id, vdata_ref, access_mode); VSgetname(vdata_id, vdata_name); nt = VFfieldtype(vdata_id, 0); if(!strcmp(name, vdata_name) && nt == DFNT_REAL){ found=TRUE; VSQuerycount(vdata_id, &n_records); VSread(vdata_id, (unsigned char*)x, n_records, FULL_INTERLACE ); VSdetach(vdata_id); break; } VSdetach(vdata_id); } } if( !found ){ printf("hget_" "real" ": ERROR searching for %s\n", name); } return found; }

int hget_float ( float *x, char *name, int vgroup_id ) { int npairs, found=FALSE, i, status; int32 vdata_tag, vdata_ref, vdata_id, n_records, nt; char vdata_name[VSNAMELENMAX]; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } npairs = Vntagrefs(vgroup_id); for(i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &vdata_tag, &vdata_ref ); if( Visvs(vgroup_id, vdata_ref) ){ vdata_id=VSattach(file_id, vdata_ref, access_mode); VSgetname(vdata_id, vdata_name); nt = VFfieldtype(vdata_id, 0); if(!strcmp(name, vdata_name) && nt == DFNT_FLOAT32){ found=TRUE; VSQuerycount(vdata_id, &n_records); VSread(vdata_id, (unsigned char*)x, n_records, FULL_INTERLACE ); VSdetach(vdata_id); break; } VSdetach(vdata_id); } } if( !found ){ printf("hget_" "float" ": ERROR searching for %s\n", name); } return found; }
int hget_double ( double *x, char *name, int vgroup_id ) { int npairs, found=FALSE, i, status; int32 vdata_tag, vdata_ref, vdata_id, n_records, nt; char vdata_name[VSNAMELENMAX]; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } npairs = Vntagrefs(vgroup_id); for(i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &vdata_tag, &vdata_ref ); if( Visvs(vgroup_id, vdata_ref) ){ vdata_id=VSattach(file_id, vdata_ref, access_mode); VSgetname(vdata_id, vdata_name); nt = VFfieldtype(vdata_id, 0); if(!strcmp(name, vdata_name) && nt == DFNT_FLOAT64){ found=TRUE; VSQuerycount(vdata_id, &n_records); VSread(vdata_id, (unsigned char*)x, n_records, FULL_INTERLACE ); VSdetach(vdata_id); break; } VSdetach(vdata_id); } } if( !found ){ printf("hget_" "double" ": ERROR searching for %s\n", name); } return found; }
# 317 "hdf_stuff.c"
int hput_int_array_1d( int_array_1d * x, char *name, int vgroup_id ){ const int32 rank=1; int32 dims[1], start[1], edges[1], base[1], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_INT32, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
int hput_float_array_1d( float_array_1d * x, char *name, int vgroup_id ){ const int32 rank=1; int32 dims[1], start[1], edges[1], base[1], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_FLOAT, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
int hput_double_array_1d( double_array_1d * x, char *name, int vgroup_id ){ const int32 rank=1; int32 dims[1], start[1], edges[1], base[1], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_DOUBLE, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }

int hput_real_array_1d( real_array_1d * x, char *name, int vgroup_id ){ const int32 rank=1; int32 dims[1], start[1], edges[1], base[1], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_REAL, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
# 365 "hdf_stuff.c"
int 
hput_int_array_2d( int_array_2d * x, char *name, int vgroup_id )
{ 
  const int32 rank=2; 
  int32 dims[2], start[2], edges[2], base[2], sds_id, istat, ref; 
  int n; 
  if (file_id == -1)
    { 
      printf("ERROR: There is no open database file!\n"); 
      return -1;
    } 
  dims[0]=x->n2; 
  dims[1]=x->n1; 
  for( n=0; n < rank; n++ ) 
    { 
      start[n]=0; 
      edges[n]=dims[n]; 
      base[n]=1;
    } 
  sds_id = SDcreate(sd_id, name, DFNT_INT32, rank, dims); 
  istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); 
  SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); 
  ref = SDidtoref(sds_id); 
  Vaddtagref(vgroup_id, DFTAG_NDG, ref ); 
  istat = SDendaccess(sds_id); 
  return TRUE;
}
int hput_float_array_2d( float_array_2d * x, char *name, int vgroup_id ){ const int32 rank=2; int32 dims[2], start[2], edges[2], base[2], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n2; dims[1]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_FLOAT, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
int hput_double_array_2d( double_array_2d * x, char *name, int vgroup_id ){ const int32 rank=2; int32 dims[2], start[2], edges[2], base[2], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n2; dims[1]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_DOUBLE, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }

int hput_real_array_2d( real_array_2d * x, char *name, int vgroup_id ){ const int32 rank=2; int32 dims[2], start[2], edges[2], base[2], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n2; dims[1]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_REAL, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
# 413 "hdf_stuff.c"
int hput_int_array_3d( int_array_3d * x, char *name, int vgroup_id ){ const int32 rank=3; int32 dims[3], start[3], edges[3], base[3], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n3; dims[1]=x->n2; dims[2]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_INT32, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
int hput_float_array_3d( float_array_3d * x, char *name, int vgroup_id ){ const int32 rank=3; int32 dims[3], start[3], edges[3], base[3], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n3; dims[1]=x->n2; dims[2]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_FLOAT, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
int hput_double_array_3d( double_array_3d * x, char *name, int vgroup_id ){ const int32 rank=3; int32 dims[3], start[3], edges[3], base[3], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n3; dims[1]=x->n2; dims[2]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_DOUBLE, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }

int hput_real_array_3d( real_array_3d * x, char *name, int vgroup_id ){ const int32 rank=3; int32 dims[3], start[3], edges[3], base[3], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n3; dims[1]=x->n2; dims[2]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_REAL, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
# 463 "hdf_stuff.c"
int hput_int_array_4d( int_array_4d * x, char *name, int vgroup_id ){ const int32 rank=4; int32 dims[4], start[4], edges[4], base[4], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n4; dims[1]=x->n3; dims[2]=x->n2; dims[3]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_INT32, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
int hput_float_array_4d( float_array_4d * x, char *name, int vgroup_id ){ const int32 rank=4; int32 dims[4], start[4], edges[4], base[4], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n4; dims[1]=x->n3; dims[2]=x->n2; dims[3]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_FLOAT, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
int hput_double_array_4d( double_array_4d * x, char *name, int vgroup_id ){ const int32 rank=4; int32 dims[4], start[4], edges[4], base[4], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n4; dims[1]=x->n3; dims[2]=x->n2; dims[3]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_DOUBLE, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }

int hput_real_array_4d( real_array_4d * x, char *name, int vgroup_id ){ const int32 rank=4; int32 dims[4], start[4], edges[4], base[4], sds_id, istat, ref; int n; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return -1; } dims[0]=x->n4; dims[1]=x->n3; dims[2]=x->n2; dims[3]=x->n1; for( n=0; n < rank; n++ ) { start[n]=0; edges[n]=dims[n]; base[n]=1; } sds_id = SDcreate(sd_id, name, DFNT_REAL, rank, dims); istat = SDwritedata(sds_id, start, NULL, edges, (unsigned char*)(x->arrayptr)); SDsetattr(sds_id, "arrayBase", DFNT_INT32, rank, base ); ref = SDidtoref(sds_id); Vaddtagref(vgroup_id, DFTAG_NDG, ref ); istat = SDendaccess(sds_id); return TRUE; }
# 522 "hdf_stuff.c"
int_array_1d * hget_int_array_1d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; int_array_1d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 1 && nt == DFNT_INT32){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_int_array_1d( dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "int" "_array_1d: ERROR searching for `%s'\n", name); } return x_; }
float_array_1d * hget_float_array_1d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; float_array_1d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 1 && nt == DFNT_FLOAT){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_float_array_1d( dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "float" "_array_1d: ERROR searching for `%s'\n", name); } return x_; }
double_array_1d * hget_double_array_1d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; double_array_1d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 1 && nt == DFNT_DOUBLE){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_double_array_1d( dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "double" "_array_1d: ERROR searching for `%s'\n", name); } return x_; }

real_array_1d * hget_real_array_1d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; real_array_1d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 1 && nt == DFNT_REAL){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_real_array_1d( dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "real" "_array_1d: ERROR searching for `%s'\n", name); } return x_; }
# 581 "hdf_stuff.c"
int_array_2d * hget_int_array_2d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; int_array_2d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 2 && nt == DFNT_INT32){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_int_array_2d( dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "int" "_array_2d: ERROR searching for `%s'\n", name); } return x_; }
float_array_2d * hget_float_array_2d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; float_array_2d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 2 && nt == DFNT_FLOAT){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_float_array_2d( dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "float" "_array_2d: ERROR searching for `%s'\n", name); } return x_; }
double_array_2d * hget_double_array_2d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; double_array_2d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 2 && nt == DFNT_DOUBLE){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_double_array_2d( dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "double" "_array_2d: ERROR searching for `%s'\n", name); } return x_; }

real_array_2d * hget_real_array_2d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; real_array_2d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 2 && nt == DFNT_REAL){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_real_array_2d( dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "real" "_array_2d: ERROR searching for `%s'\n", name); } return x_; }
# 640 "hdf_stuff.c"
int_array_3d * hget_int_array_3d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; int_array_3d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 3 && nt == DFNT_INT32 ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_int_array_3d( dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "int" "_array_3d: ERROR searching for `%s'\n", name); } return x_; }
float_array_3d * hget_float_array_3d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; float_array_3d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 3 && nt == DFNT_FLOAT ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_float_array_3d( dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "float" "_array_3d: ERROR searching for `%s'\n", name); } return x_; }
double_array_3d * hget_double_array_3d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; double_array_3d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 3 && nt == DFNT_DOUBLE ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_double_array_3d( dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "double" "_array_3d: ERROR searching for `%s'\n", name); } return x_; }

real_array_3d * hget_real_array_3d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; real_array_3d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 3 && nt == DFNT_REAL ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_real_array_3d( dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "real" "_array_3d: ERROR searching for `%s'\n", name); } return x_; }
# 699 "hdf_stuff.c"
int_array_4d * hget_int_array_4d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; int_array_4d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 4 && nt == DFNT_INT32 ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_int_array_4d( dims[3], dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "int" "_array_4d: ERROR searching for `%s'\n", name); } return x_; }
float_array_4d * hget_float_array_4d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; float_array_4d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 4 && nt == DFNT_FLOAT ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_float_array_4d( dims[3], dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "float" "_array_4d: ERROR searching for `%s'\n", name); } return x_; }
double_array_4d * hget_double_array_4d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; double_array_4d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 4 && nt == DFNT_DOUBLE ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_double_array_4d( dims[3], dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "double" "_array_4d: ERROR searching for `%s'\n", name); } return x_; }

real_array_4d * hget_real_array_4d( char * name, int vgroup_id ){ int n, i, npairs, found=FALSE, status; int32 tag, ref, rank, nt, nattrs, sds_index, sds_id , start[MAX_VAR_DIMS], edges[MAX_VAR_DIMS], dims[MAX_VAR_DIMS]; char sd_name[MAX_NC_NAME]; real_array_4d * x_ = NULL; if (file_id == -1){ printf("ERROR: There is no open database file!\n"); return NULL; } npairs = Vntagrefs(vgroup_id); for( i=0; i<npairs; i++ ){ status = Vgettagref(vgroup_id, i, &tag, &ref ); if( tag == DFTAG_NDG ){ sds_index = SDreftoindex(sd_id, ref); sds_id = SDselect(sd_id, sds_index); status = SDgetinfo(sds_id, sd_name, &rank, dims, &nt, &nattrs); if( !strcmp(name, sd_name) && rank == 4 && nt == DFNT_REAL ){ found=TRUE; for( n=0; n<rank; n++ ){ start[n]=0; edges[n]=dims[n]; } x_ = create_real_array_4d( dims[3], dims[2], dims[1], dims[0] ); status = SDreaddata(sds_id, start, NULL, edges, (unsigned char*) x_->arrayptr); status = SDendaccess(sds_id); break; } status = SDendaccess(sds_id); } } if( !found ){ printf("hget_" "real" "_array_4d: ERROR searching for `%s'\n", name); } return x_; }




int
hput_string( char * x, char * name, int vgroup_id ){
  int32 vdata_id;
  int num;
  char null_string[] = "\0";

  if (file_id == -1){
    printf("ERROR: There is no open database file!\n");
    return -1;
  }

  vdata_id = VSattach(file_id, -1, access_mode);
  VSsetname(vdata_id, name);
  VSsetclass(vdata_id, "string" );

  VSfdefine(vdata_id, name, DFNT_CHAR8, 1);

  VSsetfields(vdata_id, name);

  if (x){
    num=strlen(x)+1;
    VSwrite(vdata_id, (unsigned char*) x, num, FULL_INTERLACE);
  }
  else{
    VSwrite(vdata_id, (unsigned char*) null_string, 1, FULL_INTERLACE);
  }

  Vinsert(vgroup_id, vdata_id);

  VSdetach(vdata_id);
  return TRUE;
}


char *
hget_string(char * name, int vgroup_id ){
  int32 vdata_tag, vdata_ref, vdata_id, n_records;
  char vdata_name[VSNAMELENMAX], *x = NULL;
  int found=FALSE, npairs, status, i;

  if (file_id == -1){
    printf("ERROR: There is no open database file!\n");
    return NULL;
  }

  npairs = Vntagrefs(vgroup_id);
  for( i=0; i<npairs; i++ )
  {

    status = Vgettagref(vgroup_id, i, &vdata_tag, &vdata_ref );
    if( Visvs(vgroup_id, vdata_ref) ){

      vdata_id = VSattach(file_id, vdata_ref, access_mode);

      VSgetname(vdata_id, vdata_name);
      if( !strcmp(name, vdata_name) ){
 found=TRUE;
 VSQuerycount( vdata_id, &n_records );

        x = (char *) malloc( n_records * sizeof(char) );
 VSread(vdata_id, (unsigned char*)x, n_records, FULL_INTERLACE );
        VSdetach(vdata_id);
        break;
      }
      VSdetach(vdata_id);
    }
  }
  if( !found ){
    printf("get: ERROR searching for %s\n", name);
  }
  return x;
}
