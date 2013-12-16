#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"
#include "stdio.h"
#include "string.h"
void usage() {
	errAbort(
		"knownExonIntron - find coordinates for known exons and introns\n"
		"    Usage: knownExonIntron <refSeq-gene-list>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

static void knownExonIntron(char *geneList) {
//name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	id	name2
//NM_207014	chr1	-	67075871	67163158	67075923	67163102	10	67075871,67078739,67085754,67100417,67109640,67113051,67129424,67131499,67
//143471,67162932,	67076067,67078942,67085949,67100573,67109780,67113208,67129537,67131684,67143646,67163158,	0	WDR78
	FILE *fp;
	char buf[50000], chr[50], exonStarts[5000], exonEnds[5000], name[100],flag[50];
	int num, start, end, i;
	char *b, *e;
	fp = mustOpen(geneList, "r");
	while(fgets(buf, 50000, fp)) {
		if (buf[0] == '#')
			continue;
		sscanf(buf,"%s",flag);
		if (sscanf(buf, "%*s %s %*c %*d %*d %*d %*d %d %s %s %*s %s", chr, &num, exonStarts, exonEnds, name)!= 5)
			errAbort("error: %s", buf);
		for (i = 0, b = exonStarts, e = exonEnds; i < num; i++) {
			if (sscanf(b, "%d,", &start) != 1)
				errAbort("exon start: %s", b);
			if (sscanf(e, "%d,", &end) != 1)
				errAbort("exon end: %s", e);
			b = strchr(b, ',');
			++b;
			e = strchr(e, ',');
			++e;
			printf("@ %s %d %d %s\n", chr, start, end, name);
		}
		b = strchr(exonStarts, ',');
		++b; 
		for (i = 0, e = exonEnds; i < num - 1; i++) {
			if (sscanf(b, "%d,", &start) != 1)
				errAbort("intron end: %s", b);
			if (sscanf(e, "%d,", &end) != 1)
				errAbort("intron start: %s", e);
			b = strchr(b, ',');
			++b;
			e = strchr(e, ',');
			++e;
			printf("- %s %d %d %s\n", chr, end, start, name);
		}

	}
	fclose(fp);
} 

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	knownExonIntron(argv[1]);
	return 0;
}

