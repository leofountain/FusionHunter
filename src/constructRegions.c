#include "common.h"
#include "errabort.h"
#include "options.h"

#define D 2000	// distance between two regions
#define N 60		// minimum length of a region
#define M 10000	// maximum length of a region

void usage() {
	errAbort(
		"constructRegions - construct regions based on predicted exons.\n"
		"    Usage: constructRegions <putative_exons>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

void constructRegions(char *exons) {
	FILE *fp;
	char buf[500], chr[50], preChr[50];
	int b, e, preBeg, preEnd;

	fp = mustOpen(exons, "r");
	preChr[0] = '\0';
	preBeg = preEnd = 0;
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%s %d %d", chr, &b, &e) != 3)
			errAbort("error: %s", buf);
		if (!sameString(chr, preChr)) {
			if (preChr[0] != '\0' && preEnd - preBeg > N) 
				printf("%s %d %d\n", preChr, preBeg, preEnd);
			strcpy(preChr, chr);
			preBeg = b;
			preEnd = 0;
		}
		if (preEnd != 0 && (b - preEnd > D || preEnd - preBeg > M)) {
			if (preEnd - preBeg > N)
				printf("%s %d %d\n", chr, preBeg, preEnd);
			preBeg = b;
		}
		preEnd = e;
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	constructRegions(argv[1]);
	return 0;
}

