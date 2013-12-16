#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

struct pos {
	char chr[8];
	int beg, end;
	char orient;
};

void usage() {
	errAbort(
		"removePcrReplicate - remove paired-end reads that might be PCR duplicates.\n"
		"    Usage: removePcrReplicate <pair_output>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

static struct hash *readsHash;

//SRR018267.10000200	chr3 50271046 50271096 +	chr3 50271408 50271458 -

void removePcrReplicate(char *pairFile) {
	FILE *fp;
	char buf[500], fub[500], chr1[50], chr2[50], id[500];
	int b1, e1, b2, e2;
	
	readsHash = newHash(16);

	fp = mustOpen(pairFile, "r");
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%s %s %d %d %*c %s %d %d %*c", 
										id, chr1, &b1, &e1, chr2, &b2, &e2) != 7)
			errAbort("error: %s", buf);
		sprintf(fub, "%s:%d-%d %s:%d-%d", chr1, b1, e1, chr2, b2, e2);
		if (hashLookup(readsHash, fub))
			continue;
		hashStoreName(readsHash, fub);
		sprintf(fub, "%s:%d-%d %s:%d-%d", chr2, b2, e2, chr1, b1, e1);
		if (hashLookup(readsHash, fub))
			continue;
		hashStoreName(readsHash, fub);
		printf("%s", buf);
	}
	fclose(fp);
	hashFree(&readsHash);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	
	removePcrReplicate(argv[1]);
	
	return 0;
}

