include ../makefile.defs

SRCS= cartesian_grid.c normal_curve.c linear_interp_grid.c disc_point_mapping.c \
      hyperbolic.c generic_mapping.c tg_grid.c


OBJS= $(SRCS:.c=.o)

mappings.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	

clean:
	rm -f *.o *.a
