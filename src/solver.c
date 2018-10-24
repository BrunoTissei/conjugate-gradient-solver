/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#include "solver.h"


// cabecalho das funcoes locais para uso nesse arquivo
static void init_memory(global_t *global, value_t *value, iteration_t *iteration);

static void free_memory(global_t *global, value_t *value, iteration_t *iteration);

static void output_result(global_t *global, value_t *value, iteration_t *iteration);

static void solve(global_t *global, value_t *value, iteration_t *iteration);


/***********************
 * Executa programa principal
 ***********************/
bool run(global_t *global) {
  bool worked = TRUE;

  uint n = global->n;
  uint bandwidth = global->bandwidth;
  uint limit = bandwidth >> 1;

  value_t value;
  iteration_t iteration;
  init_memory(global, &value, &iteration);

  double **auxA = (double **) _mm_malloc((limit + 1) * sizeof(double*), 32);
  for (uint i = 0; i <= limit; ++i)
    auxA[i] = (double *) _mm_malloc(n * sizeof(double), 32);

  for (uint i = 0; i <= limit; ++i)
    worked &= generateRandomDiagonal(n, i, bandwidth, auxA[i]);

  int k = 0;
  for (uint i = 0; i < n; ++i)
    value.main[i] = auxA[0][i];

  for (uint i = 0; i < n; ++i)
    for (uint j = 1; j <= limit; ++j)
      value.A[k++] = auxA[j][i];

  for (uint i = 0; i <= limit; ++i)
    _mm_free(auxA[i]);
  _mm_free(auxA);

  // Continua com o calculo apenas se todas a diagonais foram 
  // geradas corretamente
  if (worked) {
    generateB(n, value.b);
    for (uint i = 0; i < n; ++i)
      value.x[i] = 0.0;

    solve(global, &value, &iteration);
  } else {
    fprintf(stderr, "Insert valid values.\n");
    return FALSE;
  }

  free_memory(global, &value, &iteration);
  return TRUE;
}


/***********************
 * Aloca memora necessaria para execucao do algoritmo
 ***********************/
static void init_memory(global_t *global, value_t *value, iteration_t *iteration) {
  uint n = global->n;
  uint bandwidth = global->bandwidth;
  uint max_iter = global->max_iter;

  value->main = (double *) _mm_malloc((n - 1) * 2 * sizeof(double), 32);
 
	value->A = (double *) _mm_malloc(n * bandwidth * sizeof(double), 32);
  value->b = (double *) _mm_malloc(n *           sizeof(double), 32);
  value->x = (double *) _mm_malloc(n *           sizeof(double), 32);

  iteration->norm   = (double *) _mm_malloc(max_iter * sizeof(double), 32);
  iteration->error  = (double *) _mm_malloc(max_iter * sizeof(double), 32);
  iteration->time_cg = (double *) _mm_malloc(max_iter * sizeof(double), 32);
  iteration->time_r  = (double *) _mm_malloc(max_iter * sizeof(double), 32);
}


/***********************
 * Libera toda memoria que foi alocada
 ***********************/
static void free_memory(global_t *global, value_t *value, iteration_t *iteration) {
  _mm_free(value->main);
  _mm_free(value->A);
  _mm_free(value->b);
  _mm_free(value->x);

  _mm_free(iteration->norm);
  _mm_free(iteration->error);
  _mm_free(iteration->time_cg);
  _mm_free(iteration->time_r);
  
  fclose(global->output_file);
}


/***********************
 * Imprime os resultados no arquivo especificado pelo usuario
 ***********************/
static void output_result(global_t *global, value_t *value, iteration_t *iteration) {
  uint n = global->n;
  uint max_iter = global->max_iter;

  FILE *output_file = global->output_file;

  double min_time_cg, max_time_cg, avg_time_cg;
  double min_time_r,  max_time_r,  avg_time_r;

  avg_time_cg = min_time_cg = max_time_cg = iteration->time_cg[0];
  avg_time_r  = min_time_r  = max_time_r  = iteration->time_r[0];
  
  // Calcula menor, maior e media dos tempos do metodo e de residuo
  for (uint i = 1; i < max_iter; ++i) {
    min_time_cg = min(min_time_cg, iteration->time_cg[i]);
    max_time_cg = max(max_time_cg, iteration->time_cg[i]);

    min_time_r  = min(min_time_r, iteration->time_r[i]);
    max_time_r  = max(max_time_r, iteration->time_r[i]);

    avg_time_cg += iteration->time_cg[i];
    avg_time_r  += iteration->time_r[i];
  }

  avg_time_cg /= max_iter;
  avg_time_r  /= max_iter;

  // Todos os resultados sao impressos no output_file
  fprintf(output_file, "###########\n");

  fprintf(output_file, "# CG compute time (min, avg, max): %g %g %g\n", 
      min_time_cg, avg_time_cg, max_time_cg);

  fprintf(output_file, "# Residual compute time (min, avg, max): %g %g %g\n", 
      min_time_r,  avg_time_r,  max_time_r);

  fprintf(output_file, "#\n");

  fprintf(output_file, "# Euclidian norm of residual e approximate error\n");
  for (uint i = 0; i < max_iter; ++i)
    fprintf(output_file, "# i = %d: %g %g\n", i, iteration->norm[i], iteration->error[i]);
  fprintf(output_file, "###########\n");

  fprintf(output_file, "%d\n", n);
  for (uint i = 0; i < n; ++i)
    fprintf(output_file, "%.14g ", value->x[i]);
  fprintf(output_file, "\n");
}


/***********************
 * Executa o algoritmo e chama a funcao que mostra os resultados
 ***********************/
static void solve(global_t *global, value_t *value, iteration_t *iteration) {
  conjugate_gradient(global, value, iteration);
  output_result(global, value, iteration);
}
