#ifndef input_output_h
#define input_output_h

/* data structure for i/o files */
typedef struct input_output{
  char *copy_file_name;
  FILE *read_command, *copy_command;
} input_output;

typedef struct Color {
  double r, g, b;
} color_type;

/* Marker types for X_marker */

#define DOT                 0
#define CIRCLE              1
#define BOX                 2
#define ASTERISK            3
#define CROSS               4
#define FILLED_CIRCLE       5
#define FILLED_BOX          6
#define TRIANGLE            7
#define FILLED_TRIANGLE     8
#define INV_TRIANGLE        9
#define FILLED_INV_TRIANGLE 10

#endif
