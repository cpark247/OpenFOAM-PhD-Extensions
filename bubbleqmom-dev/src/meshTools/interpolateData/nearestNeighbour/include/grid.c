#include "grid.h"
#define NDEBUG
grid *grid_create(int dims)
{
	grid *g;
	g = malloc(sizeof(grid));
	g->dims = dims;
	g->low_lim = malloc(sizeof(double) * dims);
	g->upper_lim = malloc(sizeof(double) * dims);
	g->split = malloc(sizeof(int) * dims);
	g->coords = malloc(sizeof(int) * dims);
	grid_walkReset(g);
	return g;
}

int grid_destroy(grid *g)
{
	if(g)
	{
		if(g->low_lim)
		{
			free(g->low_lim);
		}
		if(g->upper_lim)
		{
			free(g->upper_lim);
		}
		if(g->split)
		{
			free(g->split);
		}
		if(g->coords)
		{
			free(g->coords);
		}
		free(g);
		return 1;
	}
	else
	{
		return 0;
	}
}

int grid_addLim(grid *g, int dim, double lower, double upper)
{
	if(dim >= 0 && dim < g->dims)
	{
		g->low_lim[dim] = lower;
		g->upper_lim[dim] = upper;
		return 1;
	}
	else
	{
		return 0;
	}
}
void grid_print(grid *g)
{
	printf("Number of Dimensions: %i\n", g->dims);
	int i;
	for(i = 0; i < g->dims; i++)
	{
		printf("D:%i L: %lf, H: %lf, S: %i\n", i + 1, g->low_lim[i], g->upper_lim[i], g->split[i]);
	}
}
int grid_addSplit(grid *g, int dim, int split)
{
	if(dim >= 0 && dim < g->dims)
	{
		g->split[dim] = split;
		return 1;
	}
	else
	{
		return 0;
	}
}
int grid_calcPair(grid *g, int dim, int cut, double *lower, double *upper)
{
// calculate lower and upper limit for given cut number
	if(dim >= 0 && dim < g->dims)
	{
		if(cut <= g->split[dim]) // first cut number is 0
		{
			double delta = (g->upper_lim[dim] - g->low_lim[dim]) / ((double) g->split[dim] + 1.0);
			*lower = g->low_lim[dim] + (cut) * delta;
			*upper = *lower + delta;
			return 1;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return 0;
	}
}

int grid_walk_i(int pos, grid *g)
{
	if(pos < g->dims)
	{
		//check if increase is still in range
		if(g->coords[pos] + 1 <= g->split[pos])
		{
			g->coords[pos]++;
			return 1; //true
		}
		else
		{
			g->coords[pos] = 0;
			return grid_walk_i(pos + 1, g);
		}
	}
	else
	{
		return 0;    //false
	}
}
int grid_walk(grid *g)
{
	return grid_walk_i(0, g);
}
void grid_walkReset(grid *g)
{
	int i;
	for(i = 0; i < g->dims; i++)
	{
		g->coords[i] = 0;
	}
}
int grid_buildRtree(grid *g, struct Node **root, int *nBlocks)
{
	//instantiate rtree and append generated rectangles
	double low_point;
	double high_point;
	int status = 1;
	int count = 1;
	grid_walkReset(g);
	while(status)
	{
		struct Rect *trect;
		trect = malloc(sizeof(struct Rect));
		int d;
		for(d = 0; d < g->dims; d++)
		{
			grid_calcPair(g, d, g->coords[d], &low_point, &high_point);
			trect->boundary[d] = low_point;
			trect->boundary[d + g->dims] = high_point;
#ifdef DEBUG
			printf("d: %i, c: %i, Pair: low: %lf, high: %lf\n", d, g->coords[d], low_point, high_point);
#endif
		}
#ifdef DEBUG
		int i;
		int nPoints = 2* g->dims;
		for(i = 0; i < nPoints; i++)
		{
			printf("%lf,", trect->boundary[i]);
		}
		printf("\n");
#endif
		RTreeInsertRect(trect, count, root, 0);
		status = grid_walk(g);
		count++;
	}
	count--;
	*nBlocks = count;
	return 0;
}
