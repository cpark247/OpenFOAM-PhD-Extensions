#ifndef NN_H
#define NN_H
#ifdef __cplusplus
extern "C"
{
#endif
#include "grid.h"
#include "dataBlock.h"
struct nn
{
	grid *g;
	int dims;
	int solentries;
	dataBlock *db;
	struct Node *r;
	struct Rect *searchRect;
	int curBlock;
};
typedef struct nn nn;
nn *nn_create(int dim, int solentries, double *min, double *max, int *splits);
void nn_destroy(nn *n);
void nn_addData(nn *n, double *data);
void nn_getData(nn *n, double *req, double *sol);
double nn_calcSquare(double *req, double *data, int pos);
#ifdef __cplusplus
}
#endif
#endif // NN_H
