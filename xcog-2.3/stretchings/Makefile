include ../makefile.defs

SRCS= generic_stretching.c layer_stretch.c arcl_stretch.c exponential_stretch.c \
      tanh_stretch.c


OBJS= $(SRCS:.c=.o)

stretchings.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	

clean:
	rm -f *.o *.a
