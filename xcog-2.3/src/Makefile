include ../makefile.defs

SRCS= xcog.c welcome_to_xcog.c run_tutorial.c save_composite_grid.c \
      new_composite_grid.c overlap_parameters.c compute_overlap.c \
      read_curves_mappings.c save_curves_mappings.c \
      set_global_plot_mode.c set_global_boundaries.c unique_name.c \
      insert_curve.c print_mapping_list.c 


OBJS= $(SRCS:.c=.o)


../bin/$(XCOG): $(OBJS) $(UTIL)/util.a $(OVERLAP)/overlap.a $(CURVES)/curves.a $(MAPPINGS)/mappings.a $(STRETCHINGS)/stretchings.a $(TGLIB)
	if [ ! -d ../bin ]; then mkdir ../bin; fi
	$(LD) -o $@ $(CFLAGS) $(OBJS) $(CLIBS) $(XLIBS)

clean:
	rm -f *.o *.a
