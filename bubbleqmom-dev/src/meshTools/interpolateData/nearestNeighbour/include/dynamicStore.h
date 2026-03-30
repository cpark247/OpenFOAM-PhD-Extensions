#ifndef DYNAMICSTORE_H
#define DYNAMICSTORE_H

#ifndef STDLIB_H
#define STDLIB_H
#include <stdlib.h>
#endif

#ifndef STDIO_H
#define STDIO_H
#include <stdio.h>
#endif

#ifndef STRING_H
#define STRING_H
#include <string.h>
#endif

#include "mem_inlines.h"

struct dynamicStore
{
	void *data;
	size_t elementSize;
	int curNumOfElems;
	int allocatedNumOfElems;
	int defaultIncrease;
};
typedef struct dynamicStore dynamicStore;
dynamicStore *dynamicStore_create(size_t size);
void *dynamicStore_getElement(dynamicStore *ds, int num);
int dynamicStore_storeElement(dynamicStore *ds, void *d);
int dynamicStore_destroy(dynamicStore *ds);
int dynamicStore_getNumberOfElements(dynamicStore *ds);
int dynamicStore_appendToFirst(dynamicStore *dst, dynamicStore *src);
#endif
