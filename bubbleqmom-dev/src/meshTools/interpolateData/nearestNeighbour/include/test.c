#include "rtree.h"
#include <stdio.h>
#include "dataBlock.h"

struct Rect rects[] =
{
	{0, 0, 2, 2}, // xmin, ymin, xmax, ymax (for 2 dimensional RTree)
	{2, 2, 4, 4},
	{2, 0, 4, 2},
	{0, 2, 2, 4},
};
int nrects = sizeof(rects) / sizeof(rects[0]);
struct Rect search_rect =
{
	{2, 2, 4.2, 4.2}
};

int MySearchCallback(int id, void *arg)
{
	printf("Hit data rect %d\n", id);
	return 0; // keep going return 1
}

int main()
{
	struct Node *root = RTreeNewIndex();
	int i, nhits;
	printf("nrects = %d\n", nrects);
	/*
	 * Insert all the data rects.
	 * Notes about the arguments:
	 * parameter 1 is the rect being inserted,
	 * parameter 2 is its ID. NOTE: *** ID MUST NEVER BE ZERO ***, hence the +1,
	 * parameter 3 is the root of the tree. Note: its address is passed
	 * because it can change as a result of this call, therefore no other parts
	 * of this code should stash its address since it could change undernieth.
	 * parameter 4 is always zero which means to add from the root.
	 */
	for(i = 0; i < nrects; i++)
	{
		RTreeInsertRect(&rects[i], i + 1, &root, 0);    // i+1 is rect ID. Note: root can change
	}
	nhits = RTreeSearch(root, &search_rect, MySearchCallback, 0);
	printf("Search resulted in %d hits\n", nhits);
	return 0;
}
