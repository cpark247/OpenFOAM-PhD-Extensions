#include "dataBlock.h"
int main()
{
	double array[3] = {1, 2, 3};
	dataBlock *db;
	db = dataBlock_create(5, 1, sizeof(double) * 3);
	dataBlock_append(db, 2, array);
	dataBlock_print(db);
	dataBlock_destroy(db);
	return 0;
}
