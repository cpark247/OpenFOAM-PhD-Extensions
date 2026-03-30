#include "nn.h"
int main()
{
	double min[4] = {0.0, 0.0, 0.0, 0.0};
	double max[4] = {1.0, 4.0, 0.0, 1.0};
	int dims = 4;
	int splits[4] = {1, 1, 1, 1};
	nn *n = nn_create(dims, 3, min, max, splits);
	double sample[7] = {1,1.9,0,1,   1,1,1};
	double sample2[7] = {1,3,0,1,   2,2,2};
	double req[4] = {1,2.1,0,1};
	double sol[3];
	nn_addData(n, sample);
	nn_addData(n, sample2);
	nn_getData(n, req, sol);
    dataBlock_print(n->db);
	//printf("Sol: %lf\n", sol);
	nn_destroy(n);

	return 0;
}
