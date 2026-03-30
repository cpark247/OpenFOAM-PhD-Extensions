#include "nn.h"
#define NDEBUG
int MySearchCallback(int id, void *arg);
nn *nn_create(int dim, int solentries, double *min, double *max, int *splits)
{
	if(dim != NUMDIMS)
	{
		printf("Error: NUMDIMS needs to be changed in rtree.h\n");
		return NULL;
	}
	else
	{
		nn *n;
		n = malloc(sizeof(nn));
		n->dims = dim;
		n->solentries = solentries;
		n->searchRect = malloc(sizeof(struct Rect));
		n->g = grid_create(n->dims);
		int i;
		for(i = 0; i < n->dims; i++)
		{
			grid_addLim(n->g, i, min[i], max[i]);
			grid_addSplit(n->g, i, splits[i]);
		}
#ifdef DEBUG
		grid_print(n->g);
#endif
		n->r = RTreeNewIndex();
		int nBlocks;
		grid_buildRtree(n->g, &(n->r), &nBlocks);
		n->db = dataBlock_create(nBlocks, n->solentries + n->dims, sizeof(double));
		return n;
	}
}
void nn_destroy(nn *n)
{
	grid_destroy(n->g);
	dataBlock_destroy(n->db);
	free(n);
}
void nn_search_rtree(nn *n)
{
	RTreeSearch(n->r, n->searchRect, MySearchCallback, &n->curBlock);
}

int MySearchCallback(int id, void *arg)
{
#ifdef DEBUG
	printf("Found data rect %d\n", id-1);
#endif
	int *cur;
	cur = arg;
	*cur = id-1;
	return 0; // keep going return 1
}
void nn_addData(nn *n, double *data)
{
	// copy coordinates to searchRect
	int i;
	for(i = 0; i < n->dims; i++)
	{
		n->searchRect->boundary[i] = data[i];
		n->searchRect->boundary[i + n->dims] = data[i];
#ifdef DEBUG
		printf("copy to Rect: d: %i v: %lf\n", i, data[i]);
#endif
	}
	nn_search_rtree(n);
	dataBlock_append(n->db, n->curBlock, data);
#ifdef DEBUG
	printf("Found Block for adding: %i\n", n->curBlock);
	dataBlock_print(n->db);
#endif
}
void nn_getData(nn *n, double *req, double *sol)
{
	int i;
	for(i = 0; i < n->dims; i++)
	{
		n->searchRect->boundary[i] = req[i];
		n->searchRect->boundary[i + n->dims] = req[i];
#ifdef DEBUG
		printf("copy to Rect: d: %i v: %lf\n", i, req[i]);
#endif
	}
	nn_search_rtree(n);
	int lim = dataBlock_getNumberOfEntries(n->db, n->curBlock);
	double squareDist;
	double *dataRow = dataBlock_get(n->db, n->curBlock, 0);
	squareDist = nn_calcSquare(req, dataRow, n->dims);
	int smallest_pos = 0;
	for(i = 1; i < lim; i++)
	{
		dataRow = dataBlock_get(n->db, n->curBlock, i);
		double tempDist = nn_calcSquare(req, dataRow, n->dims);
		if(tempDist < squareDist)
		{
			smallest_pos = i;
			squareDist = tempDist;
		}
	}
	dataRow = dataBlock_get(n->db, n->curBlock, smallest_pos);
#ifdef DEBUG
	printf("SQDistance: %lf, pos: %i\n", squareDist, smallest_pos);
	printf("Solution: ");
#endif
	for(i = n->dims; i < n->solentries + n->dims; i++)
	{
		sol[i - n->dims] = dataRow[i];
#ifdef DEBUG
		printf("%lf,", dataRow[i]);
#endif
	}
#ifdef DEBUG
	printf("\n");
#endif
}
double nn_calcSquare(double *req, double * data, int pos)
{
	// return sum of squared coordinates till pos
	int i;
	double sum = 0.0;
	for(i = 0; i < pos; i++)
	{
		double delta = (req[i]-data[i]);
		sum+=delta*delta;
	}
#ifdef DEBUG
	printf("nn_calcSquare: %lf\n", sum);
#endif
	return sum;
}
