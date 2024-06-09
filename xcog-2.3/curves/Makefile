include ../makefile.defs

SRCS= straight_line.c smooth_poly_curve.c cubic_spline.c circular_arc.c \
      generic_curve.c

OBJS= $(SRCS:.c=.o)

curves.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	

clean:
	rm -f *.o *.a
