extern int printf(const char *format, ...);
extern int scanf(const char *format, ...);
extern size_t strlen(const char *cs);
extern int fgetc(FILE *stream);
extern int fputc(int c, FILE *stream);
extern int fflush(FILE *stream);
extern int fscanf(FILE *stream, const char *format, ...);
/*extern int sscanf( char *line, const char *format, ...);*/
extern int fclose(FILE *stream);
extern int fprintf(FILE *stream, const char *format, ...);
extern void free(void *p);
extern void *malloc(size_t size);
extern void exit(int status);
extern int strcmp( const char *cs, const char *ct );
#ifndef fileno
  extern int fileno( FILE *stream );
#endif
#ifndef LINUX
  extern int ioctl(int filedes, int request, ... );
#endif
extern void *memset( void *s, int c, size_t n );
extern int system( const char *string );
#ifndef _HPUX_SOURCE
extern time_t time( time_t *tloc );
#endif

extern double exp10(double);
#ifndef DEC_ALPHA
extern int nint(double);
#endif


