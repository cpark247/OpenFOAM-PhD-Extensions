#include "dynamicStore.h"
void *mallocAndClean(size_t size)
{
	void *p = malloc(size);
	return memset(p, 0, size);
}

dynamicStore *dynamicStore_create(size_t size)
{
	dynamicStore *ds;
	ds = malloc(sizeof(dynamicStore));
	ds->elementSize = size;
	ds->defaultIncrease = 10;
	ds->allocatedNumOfElems = ds->defaultIncrease;
	ds->data = mallocAndClean(ds->allocatedNumOfElems * ds->elementSize);
	ds->curNumOfElems = 0;
	return ds;
}

int dynamicStore_destroy(dynamicStore *ds)
{
	if(ds)
	{
		if(ds->data)
		{
			free(ds->data);
		}
		free(ds);
		return 1;
	}
	return 0;
}

void *dynamicStore_getElement(dynamicStore *ds, int num)
{
	if(num < 0)
	{
		return NULL;
	}
	if(num > ds->curNumOfElems - 1)
	{
		return NULL;
	}
	return void_ptr_add(ds->data, num * ds->elementSize);
}

int dynamicStore_storeElement(dynamicStore *ds, void *d)
{
	int status = 1;
	//check if there is enough space for current element to be stored
	if(ds->curNumOfElems >= ds->allocatedNumOfElems)
	{
		//allocate next slice and copy stuff around
		ds->allocatedNumOfElems += ds->defaultIncrease;
		void *d_new = mallocAndClean(ds->allocatedNumOfElems * ds->elementSize);
		memcpy(d_new, ds->data, ds->elementSize * ds->curNumOfElems);
		free(ds->data);
		ds->data = d_new;
	}
	void *pos = void_ptr_add(ds->data, ds->curNumOfElems * ds->elementSize);
	memcpy(pos , d, ds->elementSize);
	ds->curNumOfElems++;
	return status;
}
int dynamicStore_getNumberOfElements(dynamicStore *ds)
{
	return ds->curNumOfElems;
}
int dynamicStore_appendToFirst(dynamicStore *dst, dynamicStore *src)
{
	if(dst->elementSize != src->elementSize)
	{
		return 0;
	}
	else
	{
		int lim = dynamicStore_getNumberOfElements(src);
		int i;
		int status;
		for(i = 0; i < lim; i++)
		{
			void *val;
			val = dynamicStore_getElement(src, i);
			status = dynamicStore_storeElement(dst, val);
			if(status == 0)
			{
				return 0;
			}
		}
		return status;
	}

}
