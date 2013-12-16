#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"
#include "myutil.h"
#include "string.h"


void usage() {
	errAbort(
		"trimOriginalBwt - trim bowtie output\n"
		"    Usage: trimOriginalBwt <bwt-file>"
	);
}

static struct optionSpec options[] = {
	{"ncbi", OPTION_BOOLEAN},
	{NULL, 0},
};

static boolean ncbi;


void trimOriginalBwt(char *bwtFile) {
	FILE *fp;
	char buf[5000], fub[5000], id[500], chr[50], dna[500];
	int beg, len;
	char *ch;
	char ori;

	fp = mustOpen(bwtFile, "r");
	while (fgets(buf, 5000, fp)) {
		if (sscanf(buf, "%[^\t]\t%[^\n]", id, fub) != 2)
			errAbort("errlr: %s", buf);
		if (sscanf(fub, "%c %s %d %s %*s", &ori, chr, &beg, dna) != 4)
			errAbort("error: %s", fub);
		if (ncbi && (ch = strchr(id, ' ')))
			*ch = '\0';
		len = strlen(dna);
		if (sameString(chr, "chrM") || sameString(chr, "chrY"))
                        continue;
                if (strstr(chr,"random") != NULL || strstr(chr,"hap") != NULL || strstr(chr,"Un") != NULL)
                        continue;
		printf("%s\t%c\t%s\t%d\t%d\n", id, ori, chr, beg, len);
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	ncbi = optionExists("ncbi");
	trimOriginalBwt(argv[1]);
	return 0;
}

