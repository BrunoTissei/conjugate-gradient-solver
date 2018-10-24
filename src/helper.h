/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#ifndef _HELPER_H
#define _HELPER_H

#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

#include "solver.h"
#include "utils.h"


bool generateRandomDiagonal(uint n, uint k, uint bandwidth, double *diag);

void generateB(uint n, double *b);

double timestamp(void);

#endif
