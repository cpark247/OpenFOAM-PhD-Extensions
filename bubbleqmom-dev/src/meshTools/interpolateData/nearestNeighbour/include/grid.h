#ifndef GRID_H
#define GRID_H

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef STDIO_H
#define STDIO_H
#include <stdio.h>
#endif

#ifndef MATH_H
#define MATH_H
#include <math.h>
#endif

#include "rtree.h"

struct grid
{
	int dims;
	double *low_lim;
	double *upper_lim;
	int *split;  // number of splits in each dimension
	int *coords;  //state description
};
typedef struct grid grid;
grid *grid_create(int dims);
int grid_destroy(grid *g);
int grid_addLim(grid *g, int dim, double lower, double upper);
int grid_addSplit(grid *g, int dim, int split);
void grid_print(grid *g);
int grid_calcPair(grid *g, int dim, int cut, double *lower, double *upper);
int grid_walk_i(int pos, grid *g);
int grid_walk(grid *g);
void grid_walkReset(grid *g);
int grid_buildRtree(grid *g, struct Node **root, int *nBlocks);

#endif //GRID_H
