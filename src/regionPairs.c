#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

static struct hash *pairRegions, *regHash;

struct pos {
	char chr[8];
	int beg, end;
	char orient;
};

void usage() {
	errAbort(
		"regionPairs - find region association (i.e. candidate fusion)\n"
		"    Usage: regionPairs <regions> <pairMapping> -readOverlap=?"
	);
}

static struct optionSpec options[] = {
	{"readOverlap", OPTION_DOUBLE},
	{NULL, 0},
};

static double D;

static void createRangeTreeHash(char *fileName) {
	FILE *fp;
	char buf[500], chr[50];
	int b, e;
	struct hashEl *el;
	struct rbTree *tr;
	
	regHash = newHash(4);

	fp = mustOpen(fileName, "r");
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%s %d %d", chr, &b, &e) != 3)
			errAbort("error: %s", buf);
		if ((el = hashLookup(regHash, chr)))
			tr = (struct rbTree *)(el->val);
		else {
			tr = rangeTreeNew();
			hashAdd(regHash, chr, tr);
		}
		rangeTreeAddValCount(tr, b, e);
	}
	fclose(fp);
}

static void updateComb(int *com, char F, char R) {
// based on illumina paired-end reads
// FR						AB 
// -+	=>	0x00	-- 
// --	=>	0x01	-+
// ++	=>	0x10	+-
// +-	=>	0x11	++
	int i = 0;
	if (F == '+' && R == '+')
		i = 2; //0x10;
	if (F == '+' && R == '-')
		i = 3; //0x11;
	if (F == '-' && R == '-')
		i = 1; //0x01;
	if (F == '-' && R == '+')
		i = 0; //0x00;
	com[i] += 1;
}

static void regionPairs(char *pairMapping) {
	FILE *fp;
	char buf[500], chr1[50], chr2[50], str1[500], str2[500];
	char ori1, ori2;
	int b1, b2, e1, e2, i;
	struct hashEl *el, *lel, *rel;
	struct hashCookie cookie;
	struct rbTree *tr;
	struct range *lrg, *rrg;
	boolean lovlp, rovlp;
	int *od;

	pairRegions = newHash(8);

	fp = mustOpen(pairMapping, "r");
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%*s %s %d %d %c %s %d %d %c", chr1, &b1, &e1, &ori1, chr2, &b2, &e2, &ori2) != 8)
			errAbort("error: %s", buf);
		if (sameString(chr1, "chrM") || sameString(chr2, "chrM"))
			continue;
		lovlp = rovlp = FALSE;
		tr = (struct rbTree *)hashFindVal(regHash, chr1);
		//assert(tr);
		if (tr == NULL)
			continue;
		lrg = rangeTreeMaxOverlapping(tr, b1, e1);
		tr = (struct rbTree *)hashFindVal(regHash, chr2);
		//assert(tr);
		if (tr == NULL)
                        continue;
		rrg = rangeTreeMaxOverlapping(tr, b2, e2);
		if (lrg == NULL || rrg == NULL || lrg == rrg)
			continue;
		if (min(lrg->end, e1) - max(lrg->start, b1) + 1 > (e1 - b1) * D)
			lovlp = TRUE;
		if (min(rrg->end, e2) - max(rrg->start, b2) + 1 > (e2 - b2) * D)
			rovlp = TRUE;
		if (lovlp && rovlp) {
			sprintf(str1, "%s:%d-%d %s:%d-%d", chr1, lrg->start, lrg->end, chr2, rrg->start, rrg->end);
			sprintf(str2, "%s:%d-%d %s:%d-%d", chr2, rrg->start, rrg->end, chr1, lrg->start, lrg->end);
	
			lel = hashLookup(pairRegions, str1);
			rel = hashLookup(pairRegions, str2);
			if (lel || (lel == NULL && rel == NULL)) {
				if (lel == NULL) {
					AllocArray(od, 4);
					for (i = 0; i < 4; i++)
						od[i] = 0;
					lel = hashAdd(pairRegions, str1, od);
				}
				od = (int *)(lel->val);
				updateComb(od, ori1, ori2);
			}
			else if (rel) {
				od = (int *)(rel->val);
				updateComb(od, ori2, ori1);
			}
		}
	}
	fclose(fp);

	cookie = hashFirst(pairRegions);
	while ((el = hashNext(&cookie))) {
		printf("%s", el->name);
		for (i = 0; i < 4; i++) {
			printf("\t%d", ((int *)(el->val))[i]);
		}
		printf("\n");
	}
	hashFree(&pairRegions);
	hashFreeWithVals(&regHash, rbTreeFree);
} 

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 3)
		usage();
	optionMustExist("readOverlap");
	D = optionDouble("readOverlap", 0.5);
	createRangeTreeHash(argv[1]);
	regionPairs(argv[2]);
	return 0;
}

