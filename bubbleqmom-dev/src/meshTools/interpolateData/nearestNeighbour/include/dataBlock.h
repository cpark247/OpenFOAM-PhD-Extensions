#ifndef DATABLOCK_H
#define DATABLOCK_H
#include <stdlib.h>
#include <stdio.h>
#include "dynamicStore.h"
struct dataBlock
{
	dynamicStore **data;
	int numOfBlocks;
	int numOfElements;
	size_t elementSize;
};
typedef struct dataBlock dataBlock;
dataBlock *dataBlock_create(int num, int numOfElements, size_t elementSize);
void dataBlock_destroy(dataBlock *db);
void dataBlock_print(dataBlock *db);
int dataBlock_append(dataBlock *db, int pos, void *dataPos);
void *dataBlock_get(dataBlock *db, int pos, int element);
int dataBlock_getNumberOfEntries(dataBlock *db, int pos);
#endif
