#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

#define D 20

void usage() {
	errAbort(
		"putativeExons - find exons covered by bowtie output.\n"
		"    Usage: putativeExons <left_bwt_file> <right_bwt_file> <left_part_bwt> <right_part_bwt>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

void putativeExons(char **argv, int total) {
	FILE *fp;
	char buf[500], chr[50];
	int b, i, len;
	int *a;
	struct hash *nameHash = hashNew(8);
	struct hashEl *el;
	struct hashCookie cookie;
	struct rbTree *tr;
	struct range *rg;

	for (i = 1; i <= total; i++) {
		fp = mustOpen(argv[i], "r");
		while (fgets(buf, 500, fp)) {
			if (sscanf(buf, "%*[^\t]\t%*c %s %d %d", chr, &b, &len) != 3)	
				errAbort("error: %s", buf);
			if ((el = hashLookup(nameHash, chr)))
				tr = (struct rbTree *)(el->val);
			else {
				tr = rangeTreeNew();
				hashAdd(nameHash, chr, tr);
			}
			rangeTreeAddValCount(tr, max(0, b - D), b + len - 1 + D);
		}
		fclose(fp);
	}
	
	cookie = hashFirst(nameHash);
	while ((el = hashNext(&cookie))) {
		tr = (struct rbTree *)(el->val);
		for (rg = rangeTreeList(tr); rg; rg = rg->next) {
			a = (int *)(rg->val);
			if (*a < 5) // ignore regions covered by less than 5 reads
				continue;
			printf("%s %d %d %d\n", el->name, rg->start, rg->end, *a);
		}
	}
	hashFreeWithVals(&nameHash, rbTreeFree);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 3 && argc != 5)
		usage();
	putativeExons(argv, argc - 1);
	return 0;
}

