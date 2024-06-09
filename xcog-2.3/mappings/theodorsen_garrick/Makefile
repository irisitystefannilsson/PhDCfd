include ../../makefile.defs
CFLAGS2 = $(CFLAGS) -I../../include

.c.o:
	$(CC) -c $(CFLAGS2) $<


CFFT_SRCS = cffti.f cffti1.f cfftf.f cfftf1.f cfftb.f cfftb1.f \
            passf.f passf2.f passf3.f passf4.f passf5.f \
            passb.f passb2.f passb3.f passb4.f passb5.f
F_SRCS = tg_transform.f dcgs.f blas.f mach.f dmset.f dlaputil.f dmvops.f xersla.f 
C_SRCS= grid.c complex.c mktgcoeff.c mktggrid.c system.c

OBJS= $(C_SRCS:.c=.o) $(F_SRCS:.f=.o) $(CFFT_SRCS:.f=.o)

theodorsen_garrick.a: $(OBJS) 
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	

clean:
	rm -f *.o *.a
