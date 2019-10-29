/**
 * Copyright (c) 2018 Bruno Freitas Tissei
 *
 * Distributed under the MIT software licenser. For the full copyright and
 * license information, please view the LICENSE file distributed with this
 * source code. 
 */

#ifndef _SOLVER_H
#define _SOLVER_H

#include <stdio.h>
#include <stdlib.h>

#include "utils.h"
#include "helper.h"
#include "algorithm.h"

/**
 * Executes every procedure necessary for CG method
 *
 * @param global necessary set of variables used in cg method
 * @return whether the execution was successful or not
 */
bool run(global_t *global);

#endif
