/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * DARTS - A fine-grain dataflow-inspired runtime system.                          *
 * Copyright (C) 2011-2014  University of Delaware                                 *
 *                                                                                 *
 * This library is free software; you can redistribute it and/or                   *
 * modify it under the terms of the GNU Lesser General Public                      *
 * License as published by the Free Software Foundation; either                    *
 * version 2.1 of the License, or (at your option) any later version.              *
 *                                                                                 *
 * This library is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU               *
 * Lesser General Public License for more details.                                 *
 *                                                                                 *
 * You should have received a copy of the GNU Lesser General Public                *
 * License along with this library; if not, write to the Free Software             *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA  *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef DARTS_LULESH_MAIN_H
#define DARTS_LULESH_MAIN_H

#include "DARTS.h"
#include <stdint.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N_THREADS 1
#define TILE_SIZE 84
#define N_COLS_TILE_SZ TILE_SIZE
#define N_ROWS_TILE_SZ TILE_SIZE
#define TOTAL_TILE_SZ (N_ROWS_TILE_SZ * N_COLS_TILE_SZ)


#define N_CORES TOTAL_NUM_CU
//#define N_CORES TOTAL_NUM_CORES



#endif // DARTS_MSORT_H
