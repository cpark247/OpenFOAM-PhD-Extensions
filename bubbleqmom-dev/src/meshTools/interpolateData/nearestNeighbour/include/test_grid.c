#include "grid.h"
int main()
{
	double low[2] = {0.0, 0.0};
	double upper[2] = {1.0, 4.0};
	int splits[2] = {2, 3};
	grid *g = grid_create(2);
	int i;
	for(i = 0; i < 2; i++)
	{
		grid_addLim(g, i, low[i], upper[i]);
		grid_addSplit(g, i, splits[i]);
	}
	grid_print(g);
	grid_buildRtree(g);
	grid_destroy(g);
	return 0;
}

