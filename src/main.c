/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef LIKWID
  #include <likwid.h>
#endif

#include "helper.h"
#include "solver.h"
#include "utils.h"

/**
 * Prints correct usage of the program and optional parameters
 */
static void printUsage(void);


// Parses input and call CG method
int main(int argc, char **argv) {
  #ifdef LIKWID
    LIKWID_MARKER_INIT;
  #endif

  srand(20162);

  int n, bandwidth, max_iter;
  double tolerance;

  bool file_set = FALSE;
  FILE *output_file = NULL;

  if (argc < 3) {
    printUsage();
    return 1;
  }

  n = atoi(argv[1]);
  bandwidth = atoi(argv[2]);
  tolerance = -1.0;
  max_iter = n;

  if (bandwidth % 2 == 0) {
    fprintf(stderr, "Insert an odd number for bandwidth.\n");
    return 1;
  }

  if (n <= 0 || bandwidth <= 0) {
    fprintf(stderr, "Insert positive values for n and bandwidth.\n");
    return 1;
  }

  for (int i = 3; i < argc; i += 2) {
    switch (argv[i][1]) {
      case 'i':
        max_iter = atoi(argv[i+1]);
        break;
      case 't':
        tolerance = strtod(argv[i+1], NULL);
        break;
      case 'o':
        output_file = fopen(argv[i+1], "w+");
        file_set = TRUE;
        break;
    }
  }

  if (!file_set) {
    printUsage();
    return 0;
  }

  global_t global;

  global.n = (uint) n;
  global.bandwidth = (uint) bandwidth;
  global.max_iter = max_iter;
  global.tolerance = tolerance;
  global.output_file = output_file;

  if (!run(&global))
    return 1;

  #ifdef LIKWID
    LIKWID_MARKER_CLOSE;
  #endif

  return 0;
}


//  Prints correct usage of the program and optional parameters
static void printUsage(void) {
  fprintf(stderr, "\n----------------- cgSolver -----------------\n");
  fprintf(stderr, "usage: cgSolver <n> <bandwidth> <-o output_file> [-i max_iter] [-t tolerance]\n\n");
  fprintf(stderr, "   n: Dimension of system of linear equations (n > 0).\n");
  fprintf(stderr, "   bandwidth: Bandwidth of matrix A (must be odd).\n");
  fprintf(stderr, "   output_file: Path for output file.\n\n");

  fprintf(stderr, "   max_iter: Number of iterations (default = n).\n");
  fprintf(stderr, "   tolerance: Maximum absolute error (euclidean norm of residual).\n\n");
}
