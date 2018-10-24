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


/**
 * Generates random single diagonal, must be called (bandwidth / 2) times to 
 * build complete matrix.
 *
 * @param n linear system dimension
 * @param k number of diagonal (0: main diagonal, 1: 1st diagonal above/bellow...)
 * @param bandwidth number of bands (bandwidth) of matrix)
 * @param diag array that stores computed diagonal
 * @return whether was successful or not
 */
bool generate_random_diagonal(uint n, uint k, uint bandwidth, double *diag);


/**
 * Generates vector of independent term, where the values are given by a
 * specific generator function.
 *
 * @param n linear system dimension
 * @param b array that stores computed values
 */
void generate_b(uint n, double *b);


/**
 * Gets current time in milliseconds.
 *
 * @return current time in milliseconds
 */
double timestamp(void);

#endif
