#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

#define GAP 300
#define SZ  15

void usage() {
	errAbort(
		"chainPairs - find alignment pair in selfchains\n"
		"    Usage: chainPairs <self-chain-dir>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

void chainPairs(char *dir) {
	FILE *fp;
	char buf[500], chr[50], fileName[500], chr1[50], chr2[50];
	int id, b1, e1, gap1, gap2, i1, i2, j, ali, size, b2, e2, base1, base2;
	char ori;

	for (j = 1; j <= 23; j++) {
		if (j == 23)
			sprintf(chr, "chrX");
		else
			sprintf(chr, "chr%d", j);
		sprintf(fileName, "%s/%s.chain", dir, chr);
		fp = mustOpen(fileName, "r");
		while (fgets(buf, 500, fp)) {
			if (buf[0] == '\n')
				continue;
			if (buf[0] == 'c') {
				if (sscanf(buf, "chain %*d %s %*d %*c %d %d %s %d %c %d %d %d",
									 chr1, &b1, &e1, chr2, &size, &ori, &b2, &e2, &id) != 9) 
					errAbort("chain: %s", buf);
				if (ori == '+') {
					i1 = base1 = b1;
					i2 = base2 = b2;
					for (;;) {
						if (!fgets(buf, 500, fp))
							errAbort("%s error", fileName);
						if (sscanf(buf, "%d %d %d", &ali, &gap1, &gap2) == 3) {
							if (gap1 > GAP || gap2 > GAP) {
								if (i1 + ali - base1 > SZ) {
									printf("%s %d %d %s %d %d + [%d]\n", chr1, base1 + 1,  i1 + ali, chr2, base2 + 1, i2 + ali, id);
									printf("%s %d %d %s %d %d + [%d]\n", chr2, base2 + 1,  i2 + ali, chr1, base1 + 1, i1 + ali, id);
								}
								base1 = i1 + ali + gap1;
								base2 = i2 + ali + gap2; 
							}
							i1 += (ali + gap1);
							i2 += (ali + gap2);
						}
						else if (sscanf(buf, "%d", &ali) == 1) {
							if (i1 + ali - base1 > SZ) {
								printf("%s %d %d %s %d %d + [%d]\n", chr1, base1 + 1, i1 + ali, chr2, base2 + 1, i2 + ali, id);
								printf("%s %d %d %s %d %d + [%d]\n", chr2, base2 + 1, i2 + ali, chr1, base1 + 1, i1 + ali, id);
							}
							break;
						}
						else errAbort("chain: %s", buf);
					}
				}
				else {
					i1 = base1 = b1;
					i2 = base2 = b2;
					for (;;) {
						if (!fgets(buf, 500, fp))
							errAbort("%s error", fileName);
						if (sscanf(buf, "%d %d %d", &ali, &gap1, &gap2) == 3) {
							if (gap1 > GAP || gap2 > GAP) {
								if (i1 + ali - base1 > SZ) {
									printf("%s %d %d %s %d %d - [%d]\n", chr1, base1 + 1, i1 + ali, chr2, size - (i2 + ali) + 1, size - base2, id);
									printf("%s %d %d %s %d %d - [%d]\n", chr2, size - (i2 + ali) + 1, size - base2, chr1, base1 + 1, i1 + ali, id);
								}
								base1 = i1 + ali + gap1;
								base2 = i2 + ali + gap2;
							}
							i1 += (ali + gap1);
							i2 += (ali + gap2);
						}
						else if (sscanf(buf, "%d", &ali) == 1) {
							if (i1 + ali - base1 > SZ) {
								printf("%s %d %d %s %d %d - [%d]\n", chr1, base1 + 1, i1 + ali, chr2, size - (i2 + ali) + 1, size - base2, id);
								printf("%s %d %d %s %d %d - [%d]\n", chr2, size - (i2 + ali) + 1, size - base2, chr1, base1 + 1, i1 + ali, id);
							}
							break;
						}
						else errAbort("chain: %s", buf);
					}
				}
			}
		}
		fclose(fp);
	}
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	chainPairs(argv[1]);
	return 0;
}

