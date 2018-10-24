/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#include "helper.h"


// Generates random single diagonal, must be called (bandwidth / 2) times to 
// build complete matrix.
bool generate_random_diagonal(uint n, uint k, uint bandwidth, double *diag) {
  if (!diag || n < 3 || bandwidth > n / 2 || k > bandwidth)
    return FALSE;

  double factor = (k == 0) ? ((double) (bandwidth - 1)) : (0.0);
  double inv_rand_max = 1.0 / (double) RAND_MAX;

  for (uint i = 0; i < n - k; ++i)
    diag[i] = factor + (double) rand() * inv_rand_max;

  return TRUE;
}


// Generates vector of independent term, where the values are given by a
// specific generator function.
void generate_b(uint n, double *b) {
  double x;
  double pi = M_PI;
  double pi_x_2 = 2 * pi;
  double pi2_x_4 = 2 * pi_x_2 * pi;
  double accum_pi = 0.0;

  for (uint i = 0; i < n; ++i) {
    x = accum_pi / ((double) n);
    b[i] = pi2_x_4 * (sin(pi_x_2 * x) + sin(pi_x_2 * (pi - x)));
    accum_pi += pi;
  }
}


// Gets current time in milliseconds.
double timestamp(void) {
  struct timeval tp;
  gettimeofday(&tp, NULL);

  return((double) (tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}
