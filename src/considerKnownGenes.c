#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

#define D 50 // extention to both sides of a region

void usage() {
	errAbort(
		"considerKnownGenes - make sure that an annotated known gene is not in two separate regions.\n"
		"    Usage: considerKnownGenes <original_region> <known_gene>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

void considerKnownGenes(char **argv) {
	FILE *fp;
	char buf[10000], chr[50];
	int b, e, i;
	struct hash *nameHash = hashNew(8);
	struct hashEl *el;
	struct hashCookie cookie;
	struct rbTree *tr;
	struct range *rg;

	for (i = 1; i <= 2; i++) {
		fp = mustOpen(argv[i], "r");
		while (fgets(buf, 10000, fp)) {
			if (buf[0] == '#')
				continue;
			if (i == 1) {
				if (sscanf(buf, "%s %d %d", chr, &b, &e) != 3)	
					errAbort("error: %s", buf);
			}
			if (i == 2) {
				if (sscanf(buf, "%*s %s %*c %d %d %*s", chr, &b, &e) != 3)	
					errAbort("error: %s", buf);
			}
			if ((el = hashLookup(nameHash, chr)))
				tr = (struct rbTree *)(el->val);
			else {
				tr = rangeTreeNew();
				hashAdd(nameHash, chr, tr);
			}
			if (i == 1)
				rangeTreeAddValCount(tr, max(0, b - D), e + D);
			if (i == 2 && rangeTreeOverlaps(tr, b, e)) 
				rangeTreeAddValCount(tr, max(0, b - D), e + D);
		}
		fclose(fp);
	}
	
	cookie = hashFirst(nameHash);
	while ((el = hashNext(&cookie))) {
		tr = (struct rbTree *)(el->val);
		for (rg = rangeTreeList(tr); rg; rg = rg->next) {
			printf("%s %d %d\n", el->name, rg->start, rg->end);
		}
	}
	hashFreeWithVals(&nameHash, rbTreeFree);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 3)
		usage();
	considerKnownGenes(argv);
	return 0;
}

