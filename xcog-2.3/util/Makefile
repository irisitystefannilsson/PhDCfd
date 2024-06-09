include ../makefile.defs

SRCS= c_array.c spline.c nrutil.c xplot.c min_max.c general_help.c print_help.c \
      open_file.c log_cosh.c hdf_stuff.c linked_list.c

OBJS= $(SRCS:.c=.o)

util.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@	

clean:
	rm -f *.o *.a
