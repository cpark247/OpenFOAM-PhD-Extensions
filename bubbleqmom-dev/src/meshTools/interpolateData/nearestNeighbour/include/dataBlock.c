#include "dataBlock.h"
dataBlock *dataBlock_create(int num, int numOfElements, size_t elementSize)
{
	dataBlock *db = malloc(sizeof(dataBlock));
	db->numOfBlocks = num;
	db->numOfElements = numOfElements;
	db->elementSize = elementSize;
	db->data = malloc(sizeof(dynamicStore *) * db->numOfBlocks);
	int i;
	size_t datasize = db->numOfElements * db->elementSize;
	for(i = 0; i < db->numOfBlocks; i++)
	{
		db->data[i] = dynamicStore_create(datasize);
	}
	return db;
}
void dataBlock_destroy(dataBlock *db)
{
	if(db)
	{
		if(db->data)
		{
			int i;
			for(i = 0; i < db->numOfBlocks; i++)
			{
				dynamicStore_destroy(db->data[i]);
			}
			free(db->data);
		}
		free(db);
	}
}
void dataBlock_print(dataBlock *db)
{
	int i;
	for(i = 0; i < db->numOfBlocks; i++)
	{
		int lim = dynamicStore_getNumberOfElements(db->data[i]);
		printf("Block: %i\t", i);
		printf("Elements: %i\n", lim);
	}
}
int dataBlock_append(dataBlock *db, int pos, void *dataPos)
{
	return dynamicStore_storeElement(db->data[pos], dataPos);
}
void *dataBlock_get(dataBlock *db, int pos, int element)
{
	return  dynamicStore_getElement(db->data[pos], element);
}
int dataBlock_getNumberOfEntries(dataBlock *db, int pos)
{
	return dynamicStore_getNumberOfElements(db->data[pos]);
}

