/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#include "helper.h"


/***********************
 * Gera diagonal aleatoria, deve ser chamado bandwidth / 2 vezes para representar a matriz
 *
 * N: tamanho do sistema linear
 * k: numero da diagonal, 0 = diagonal principal, 1 = acima/abaixo da diagonal, 2 = ...
 * kMax: numero de bandas do sistema linear
 * diag: vetor para armazenar os valores da diagonal. Deve ser alocado por quem chama a função.
 ***********************/
bool generateRandomDiagonal(uint n, uint k, uint bandwidth, double *diag) {
  if (!diag || n < 3 || bandwidth > n / 2 || k > bandwidth)
    return FALSE;

  // garante valor dominante para diagonal principal
  double factor = (k == 0) ? ((double) (bandwidth - 1)) : (0.0);
  double inv_rand_max = 1.0 / (double) RAND_MAX;

  for (uint i = 0; i < n - k; ++i)
    diag[i] = factor + (double) rand() * inv_rand_max;

  return TRUE;
}


/***********************
 * Gera vetor de termos independentes com a funcao fornecida
 *
 * N: tamanho do sistema linear
 * b: vetor de termos independentes do sistema linear  
 ***********************/
void generateB(uint n, double *b) {
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


/***********************
 * Retorna tempo em milisegundos atual
 ***********************/
double timestamp(void) {
  struct timeval tp;
  gettimeofday(&tp, NULL);

  return((double) (tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}
