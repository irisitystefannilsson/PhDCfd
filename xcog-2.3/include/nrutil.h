void 
nrerror(char error_text[]);
real *
vector(long nl, long nh);
void 
free_vector(real *v, long nl, long nh);
void 
read_vector( int fd, real *x, int low, int high );
void 
write_vector( int fd, real *x, int low, int high );
int 
hput_vector( real x[], int low, int high, char *name, int vgroup_id );
real * 
hget_vector( char * name, int *low, int *high, int vgroup_id );



