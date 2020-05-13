#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


void print_help() {
  fprintf(stderr, "Usage: \n  -i [input file path] -o [output file path]  -l )\n");
  exit(EXIT_FAILURE);
}


int main(int argc, char *argv[])
{
  bool do_log = false;
  char *input = NULL;
  char *output= NULL;

  if (argc < 3)
    print_help();

  for (size_t n = 1; n < argc; n++) {
    if (argv[n][0] == '-') {
      switch (argv[n][1]) {
        case 'i': input = argv[n + 1]; break;
        case 'o': output = argv[n + 1]; break;
        case 'l': do_log = true; break;
        case 'h':
        default:
          print_help();
      }
    }
  }

  FILE *fi = fopen(input, "r");
  FILE *fo = fopen(output, "w");
  char *line = NULL;
  long row_count = 0;
  size_t buf_len = 0;
  ssize_t line_len;

  if (fi == NULL || fo == NULL)
    exit(EXIT_FAILURE);

  getline(&line, &buf_len, fi);
  fputs(do_log ? "log_length\n" : "length\n", fo);
  while ((line_len = getline(&line, &buf_len, fi)) != -1) {
    if (do_log)
      fprintf(fo, "%f\n", log((double)line_len));
    else
      fprintf(fo, "%zu\n", line_len);
    if (++row_count % 1000 == 0)
      printf("%zu rows\n", row_count);
  }

  printf("Done\n");
  fclose(fi);
  fclose(fo);

  if (line)
    free(line);

  exit(EXIT_SUCCESS);
}
