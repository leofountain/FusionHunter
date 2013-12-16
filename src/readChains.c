#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

void usage() {
	errAbort(
		"readChains - find genomic segments that overlap with selfchains\n"
		"    Usage: readChains <self-chain-dir>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

void readChains(char *dir) {
	FILE *fp;
	char buf[500], chr[50], fileName[500], ref[50];
	int b, e, gap, i, j, ali;
	struct rbTree *tr;
	struct range *rg;

	for (j = 1; j <= 23; j++) {
		if (j == 23)
			sprintf(chr, "chrX");
		else
			sprintf(chr, "chr%d", j);
		sprintf(fileName, "%s/%s.chain", dir, chr);
		tr = rangeTreeNew();
		fp = mustOpen(fileName, "r");
		while (fgets(buf, 500, fp)) {
			if (buf[0] == '\n')
				continue;
			if (buf[0] == 'c') {
				if (sscanf(buf, "chain %*d %s %*d %*c %d %d %*s", ref, &b, &e) != 3) 
					errAbort("chain: %s", buf);
				assert(sameString(ref, chr));
				i = b;
				for (;;) {
					if (!fgets(buf, 500, fp))
						errAbort("%s error", fileName);
					if (sscanf(buf, "%d %d %*d", &ali, &gap) == 2) {
						rangeTreeAdd(tr, i, i + ali);
						i += (ali + gap);
					}
					else if (sscanf(buf, "%d", &ali) == 1) {
						rangeTreeAdd(tr, i, i + ali);
						assert(i + ali == e);
						break;
					}
					else errAbort("chain: %s", buf);
				}
			}
		}
		fclose(fp);
		for (rg = rangeTreeList(tr); rg; rg = rg->next) 
			printf("%s %d %d\n", chr, rg->start, rg->end);
		rbTreeFree(&tr);
	}
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	readChains(argv[1]);
	return 0;
}

