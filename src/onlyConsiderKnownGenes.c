#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

#define D 50 // extention to both sides of a region

void usage() {
	errAbort(
		"onlyConsiderKnownGenes - only retain regions that contain at least one known genes.\n"
		"    Usage: onlyConsiderKnownGenes <original_region> <known_gene>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

void onlyConsiderKnownGenes(char *regions, char *geneList) {
	FILE *fp;
	char buf[10000], chr[50];
	int b, e;
	struct hash *nameHash = hashNew(8);
	struct hashEl *el;
	struct rbTree *tr;

	fp = mustOpen(geneList, "r");
	while (fgets(buf, 10000, fp)) {
		if (buf[0] == '#')
			continue;
		if (sscanf(buf, "%*s %s %*c %d %d %*s", chr, &b, &e) != 3)	
			errAbort("error: %s", buf);
		if ((el = hashLookup(nameHash, chr)))
			tr = (struct rbTree *)(el->val);
		else {
			tr = rangeTreeNew();
			hashAdd(nameHash, chr, tr);
		}
		rangeTreeAddValCount(tr, b, e);
	}
	fclose(fp);

	fp = mustOpen(regions, "r");
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%s %d %d", chr, &b, &e) != 3)	
			errAbort("error: %s", buf);
		tr = (struct rbTree *)hashMustFindVal(nameHash, chr);
		if (rangeTreeOverlaps(tr, b, e))
			printf("%s %d %d\n", chr, b, e);
	}
	fclose(fp);

	hashFreeWithVals(&nameHash, rbTreeFree);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 3)
		usage();
	onlyConsiderKnownGenes(argv[1], argv[2]);
	return 0;
}

