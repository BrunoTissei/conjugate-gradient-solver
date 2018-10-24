/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#include "algorithm.h"


/**
 * Computes multiplication of band matrix (A) and vector (x) (Ax = y).
 *
 * @param global necessary set of variables used in cg method (n and bandwidth)
 * @param value arrays and matrix composing linear system (main diagonal, A)
 * @param v vector multiplied by matrix 
 * @param ans resulting vector
 */
static inline void multiply_matrix_vector(global_t *global, value_t *value, double *v, double *ans) {
  #ifdef LIKWID
    LIKWID_MARKER_START("MMV");
  #endif

  uint i, j, k, c;
  uint n = global->n;
  uint limit = (global->bandwidth >> 1);

  __m256d tmat, tmul;
  __m256d tvec1, tvec2, tans1, tans2;

  // Multiply main diagonal
  for (i = 0; i < n - (n % 4); i += 4) {
    tmat = _mm256_load_pd(value->main + i);
    tvec1 = _mm256_load_pd(v + i);
    tmul = _mm256_mul_pd(tmat, tvec1);

    _mm256_stream_pd(ans + i, tmul);
  }

  for ( ; i < n; ++i)
    ans[i] = value->main[i] * v[i];


  // Multiply rest of matrix
  for (i = k = 0, c = 1; i < n - 1; ++i, c -= (limit - 1)) {
    tans1 = _mm256_setzero_pd();
    tvec2 = _mm256_set_pd(v[i], v[i], v[i], v[i]);

    for (j = 0; j < limit - (limit % 4); j += 4, c += 4, k += 4) {
      tmat = _mm256_loadu_pd(value->A + k);

      tvec1 = _mm256_load_pd(v + c);
      tans1 = _mm256_fmadd_pd(tmat, tvec1, tans1);
      tans2 = _mm256_load_pd(ans + c);
      tans2 = _mm256_fmadd_pd(tmat, tvec2, tans2);

      _mm256_storeu_pd(ans + c, tans2);
    }

    tans1 = _mm256_hadd_pd(tans1, tans1);
    ans[i] += ((double*) &tans1)[0] + ((double*) &tans1)[2];

    for (; j < limit; ++j, ++c, ++k) {
      ans[i] += value->A[k] * v[c];
      ans[c] += value->A[k] * v[i];
    }
  }

  #ifdef LIKWID
    LIKWID_MARKER_STOP("MMV");
  #endif
}


/**
 * Computes dot product between two vectors.
 *
 * @param n size of vectors
 * @param a vector a
 * @param b vector b
 * @return dot product
 */
static inline double dot_product(uint n, double *a, double *b) {
  #ifdef LIKWID
    LIKWID_MARKER_START("MVV");
  #endif

  uint i;
  double ans;

  __m256d ta, tb, res;
  __m256d sum = _mm256_setzero_pd();

  for (i = 0; i < n - (n % 4); i += 4) {
    ta = _mm256_load_pd(a + i);
    tb = _mm256_load_pd(b + i);
    sum = _mm256_fmadd_pd(ta, tb, sum);
  }

  res = _mm256_hadd_pd(sum, sum);
  ans = ((double*) &res)[0] + ((double*) &res)[2];

  for ( ; i < n; ++i)
    ans += a[i] * b[i];

  #ifdef LIKWID
    LIKWID_MARKER_STOP("MVV");
  #endif

  return ans;
}


/**
 * Computes subtraction between two vectors.
 *
 * @param n size of vectors
 * @param ans resulting vector (a - b)
 * @param a vector a
 * @param b vector b
 */
static inline void subtract_vector(uint n, double *ans, double *a, double *b) {
  uint i;

  __m256d ta, tb, res;

  for (i = 0 ; i < n - (n % 4); i += 4) {
    ta = _mm256_load_pd(a + i);
    tb = _mm256_load_pd(b + i);
    res = _mm256_sub_pd(ta, tb);

    _mm256_stream_pd(ans + i, res);
  }

  for ( ; i < n; ++i)
    ans[i] = a[i] - b[i];
}


/**
 * Computes addition between a vector and a scaled one.
 *
 * @param n size of vectors
 * @param ans resulting vector (a + b*scalar)
 * @param a vector a
 * @param b vector b to be scaled
 * @param scalar number that will scale vector b
 */
static inline void add_scaled_vector(uint n, double *ans, double *a, double *b, double scalar) {
  uint i;

  __m256d sca = _mm256_set_pd(scalar, scalar, scalar, scalar);
  __m256d tans, ta, tb;

  for (i = 0; i < n - (n % 4); i += 4) {
    ta = _mm256_load_pd(a + i);
    tb = _mm256_load_pd(b + i);
    tans = _mm256_fmadd_pd(tb, sca, ta);

    _mm256_stream_pd(ans + i, tans);
  }

  for ( ; i < n; ++i)
    ans[i] = a[i] + (b[i] * scalar);
}


/**
 * Attributes a vector to another one.
 *
 * @param n size of vectors
 * @param vec vector that receives values
 * @param val vector containing values
 */
static inline void attribute_vector(uint n, double *vec, double *val) {
  uint i;

  __m256d tval;

  for (i = 0; i < n - (n % 4); i += 4) {
    tval = _mm256_load_pd(val + i);
    _mm256_stream_pd(vec + i, tval);
  }

  for ( ; i < n; ++i)
    vec[i] = val[i];
}


// Computes solution of linear system (Ax = b) using the conjugate gradient method
// and fills value->x with answer.
void conjugate_gradient(global_t *global, value_t *value, iteration_t *iteration) {
  uint i = 0;
  uint n = global->n;
  uint max_iter = global->max_iter;

  double *b = value->b;
  double *x = value->x;
  double tolerance = global->tolerance;

  // Necessary arrays for method
  double *r = (double *) _mm_malloc(n * sizeof(double), 32);
  double *v = (double *) _mm_malloc(n * sizeof(double), 32);
  double *z = (double *) _mm_malloc(n * sizeof(double), 32);

  double curr_iter, past_iter;
  double alfa, error;
  double time_cg, time_r;

  attribute_vector(n, r, b);
  attribute_vector(n, v, r);

  // Error obtained from euclidean norm of residual
  curr_iter = dot_product(n, r, r);
  error = sqrt(curr_iter);

  for ( ; i < max_iter && error > tolerance; ++i) {
    time_cg = timestamp();

    #ifdef LIKWID
      LIKWID_MARKER_START("CG");
    #endif

    multiply_matrix_vector(global, value, v, z);
    alfa = curr_iter / dot_product(n, v, z);
    add_scaled_vector(n, x, x, v, alfa);

    time_r = timestamp();

    // Either compute exact residual
    if (curr_iter < 1e-60) {
      multiply_matrix_vector(global, value, x, r);
      subtract_vector(n, r, b, r);

    // Or adjust what was already computed
    } else
      add_scaled_vector(n, r, r, z, -alfa);

    // Store compute time of residual
    iteration->time_r[i] = timestamp() - time_r;

    past_iter = curr_iter;
    curr_iter = dot_product(n, r, r);

    add_scaled_vector(n, v, r, v, (curr_iter / past_iter));

    // Computes relative approximate error
    error = fabs(sqrt(curr_iter) - sqrt(past_iter));

    #ifdef LIKWID
      LIKWID_MARKER_STOP("CG");
    #endif

    // Store compute of iteration
    iteration->time_cg[i] = (timestamp() - time_cg);

    // Store error and residual norm per iteration
    iteration->norm[i]  = sqrt(curr_iter);
    iteration->error[i] = error;
  }

  global->max_iter = i;

  _mm_free(r);
  _mm_free(v);
  _mm_free(z);
}
