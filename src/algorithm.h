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

void conjugate_gradient(global_t *global, value_t *value, iteration_t *iteration);

#endif
