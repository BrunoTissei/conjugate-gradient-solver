/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#ifndef _ALGORITHM_H
#define _ALGORITHM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>

#ifdef LIKWID
  #include <likwid.h>
#endif

#include "utils.h"
#include "helper.h"


#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)


/**
 * Computes solution of linear system (Ax = b) using the conjugate gradient method
 * and fills value->x with answer.
 *
 * @param global necessary set of variables used in cg method
 * @param value arrays and matrix composing linear system
 * @param iteration arrays of metadata to describe iterations
 */
void conjugate_gradient(global_t *global, value_t *value, iteration_t *iteration);

#endif
