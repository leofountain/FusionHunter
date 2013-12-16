#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"
#include "myutil.h"

void usage() {
	errAbort(
		"reduceBwt - eliminate spurious reads\n"
		"    Usage: reduceBwt <bwt-file> <repeats> <chain-pairs> -repeatOverlap=? -chainNum=?"
	);
}

static struct optionSpec options[] = {
	{"repeatOverlap", OPTION_DOUBLE},
	{"chainNum", OPTION_INT},
	{"ncbi", OPTION_BOOLEAN},
	{NULL, 0},
};

struct hash *chainHash = NULL, *repeatHash = NULL;
static int NUM;
static double D;
static boolean ncbi;

static void readChainPairs(char *chainFile) {
	FILE *fp;
	char buf[500], chr1[50], chr2[50];
	int b1, e1, b2, e2; 
	struct hashEl *el;
	struct rbTree *tr;

	fp = mustOpen(chainFile, "r");
	chainHash = newHash(8);
	while(fgets(buf, 500, fp)) {
		if (sscanf(buf, "%s %d %d %s %d %d", chr1, &b1, &e1, chr2, &b2, &e2) != 6)
			errAbort("error: %s", buf);
		
		el = (hashLookup(chainHash, chr1));
		if (el == NULL) {
			tr = rangeTreeNew();
			hashAdd(chainHash, chr1, tr);
		}
		else
			tr = (struct rbTree *)(el->val);
		rangeTreeAddValCount(tr, b1, e1);
	
		el = (hashLookup(chainHash, chr2));
		if (el == NULL) {
			tr = rangeTreeNew();
			hashAdd(chainHash, chr2, tr);
		}
		else
			tr = (struct rbTree *)(el->val);
		rangeTreeAddValCount(tr, b2, e2);

		if(!fgets(buf, 500, fp))
			errAbort("error: %s", chainFile);
	}
	fclose(fp);
}

static void readRepeats(char *repeatFile) {
	FILE *fp;
	char buf[500], chr[50];
	int b, e; 
	struct hashEl *el;
	struct rbTree *tr;

	fp = mustOpen(repeatFile, "r");
	repeatHash = newHash(8);
	while(fgets(buf, 500, fp)) {
		if (buf[0] == '#')
			continue;
		if (sscanf(buf, "%s %d %d %*s", chr, &b, &e) != 3)
			errAbort("error: %s", buf);
		el = (hashLookup(repeatHash, chr));
		if (el == NULL) {
			tr = rangeTreeNew();
			hashAdd(repeatHash, chr, tr);
		}
		else
			tr = (struct rbTree *)(el->val);
		rangeTreeAdd(tr, b, e);
	}
	fclose(fp);
}

void reduceBwt(char *bwtFile, char *repeatFile, char *chainFile) {
	FILE *fp;
	char buf[5000], fub[5000], id[500], chr[50];
	int beg, len;
	int *a;
	char *ch;
	struct rbTree *tr;
	struct range *rg;

	readRepeats(repeatFile);
	readChainPairs(chainFile);

	fp = mustOpen(bwtFile, "r");
	while (fgets(buf, 5000, fp)) {
		if (sscanf(buf, "%[^\t]\t%[^\n]", id, fub) != 2)
			errAbort("errlr: %s", buf);
		if (sscanf(fub, "%*c %s %d %d", chr, &beg, &len) != 3)
			errAbort("error: %s", fub);
		if (sameString(chr, "chrM") || sameString(chr, "chrY"))
			continue;
		if (strstr(chr,"random") != NULL || strstr(chr,"hap") != NULL || strstr(chr,"Un") != NULL)
			continue;
		if (ncbi && (ch = strchr(id, ' ')))
			*ch = '\0';
		tr = (struct rbTree *)hashFindVal(repeatHash, chr);
		if(tr != NULL){
		rg = rangeTreeMaxOverlapping(tr, beg, beg + len);
		if (rg && overlap(rg->start, rg->end, beg, beg + len) > D * len)
			continue;
		}
	
		tr = (struct rbTree *)hashFindVal(chainHash, chr);
		if(tr != NULL){
		rg = rangeTreeMaxOverlapping(tr, beg, beg + len);
		if (rg) {
			a = (int *)(rg->val);
			if (*a > NUM && overlap(rg->start, rg->end, beg, beg + len) > D * len)
				continue;
		}
		}
		printf("%s\t%s\n", id, fub);
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);

	//if (argc != 4)
	//	usage();
	optionMustExist("repeatOverlap");
	optionMustExist("chainNum");
	NUM = optionInt("chainNum", 10);
	D = optionDouble("repeatOverlap", 0.5);
	ncbi = optionExists("ncbi");
	reduceBwt(argv[1], argv[2], argv[3]);

	hashFreeWithVals(&chainHash, rbTreeFree);
	hashFreeWithVals(&repeatHash, rbTreeFree);
	return 0;
}

