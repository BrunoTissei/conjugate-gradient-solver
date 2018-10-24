/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#ifndef _UTILS_H
#define _UTILS_H

#define TRUE 1
#define FALSE 0

typedef unsigned int uint;
typedef unsigned short bool;


typedef struct value_t {
  double *main;
  double *A;
  double *b;
  double *x;
} value_t;


typedef struct global_t {
  double tolerance;
  FILE *output_file;
  uint n, bandwidth, max_iter;
} global_t;


typedef struct iteration_t {
  double *norm;
  double *error;
  double *time_cg;
  double *time_r;
} iteration_t;

#endif
