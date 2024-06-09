include ../makefile.defs

SRCS= overlap.c overlapping_grid.c bi_linear.c bi_cubic.c

OBJS= $(SRCS:.c=.o)

overlap.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	

clean:
	rm -f *.o *.a
