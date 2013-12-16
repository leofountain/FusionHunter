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
		"leftRightOvlp - retain paired-end reads whose left and right can be uniquely mapped\n"
		"    Usage: leftRightOvlp <left_bwt_output> <right_bwt_out> <repeats> [-ncbi] -repeatOverlap=?"
	);
}

static struct optionSpec options[] = {
	{"ncbi", OPTION_BOOLEAN},
	{"repeatOverlap", OPTION_DOUBLE},
	{NULL, 0},
};

static struct hash *leftReads, *homoHash;
static boolean ncbi = FALSE;
static double D;

static void createRangeTreeHash(char *fileName) {
	FILE *fp;
	char buf[500], chr[50];
	int b, e;
	struct hashEl *el;
	struct rbTree *tr;
	
	homoHash = newHash(8);

	fp = mustOpen(fileName, "r");
	while (fgets(buf, 500, fp)) {
		if (buf[0] == '#')
			continue;
		if (sscanf(buf, "%s %d %d", chr, &b, &e) != 3)
			errAbort("error: %s", buf);
		if ((el = hashLookup(homoHash, chr)))
			tr = (struct rbTree *)(el->val);
		else {
			tr = rangeTreeNew();
			hashAdd(homoHash, chr, tr);
		}
		rangeTreeAddValCount(tr, b, e);
	}
	fclose(fp);
}

static void createReadsHash(char *fileName) {
	FILE *fp;
	char buf[500], chrom[20], id[500];
	char ori;
	char *ch;
	int b, len;
	struct pos *pp;
	
	leftReads = newHash(8);
	fp = mustOpen(fileName, "r");
	while (fgets(buf, 500, fp)) {
		//	if (ncbi) {
		if (sscanf(buf, "%[^\t]\t%c %s %d %d", id, &ori, chrom, &b, &len) != 5)
			errAbort("error: %s", buf);
		if ((ch = strchr(id, ' ')))
			*ch = '\0';
		//		printf("%s\n",id);
		//	}
		//	else {
		//	if (sscanf(buf, "%[^\t]\t%c %s %d %d", id, &ori, chrom, &b, &len) != 5)
		//			errAbort("error: %s", buf);
		//	}
		AllocVar(pp);
		strcpy(pp->chr, chrom);
		pp->orient = ori;
		pp->beg = b;
		pp->end = b + len - 1;
//		printf("%d\t%c\t%d\n", b, ori, len);
		hashAdd(leftReads, id, pp);
		//		printf("kkkkkkkkkk\n");
	//	printf("%d\t%c\t%d\n", b, ori, len);
	}
	fclose(fp);
}

void leftRightOvlp(char *rightReads) {
	FILE *fp;
	char buf[500], fub[500], chrom[20], id[500];
	char ori;
	char *ch;
	struct pos *pp;
	struct hashEl *el;
	struct rbTree *tr;
	struct range *lrg, *rrg;
	int len, b, l;
	boolean lovlp, rovlp;
	fp = mustOpen(rightReads, "r");
	while (fgets(buf, 500, fp)) {
	//	if (ncbi) {
			if (sscanf(buf, "%[^\t]\t%c %s %d %d", id, &ori, chrom, &b, &l) != 5)
				errAbort("error: %s", buf);
			if ((ch = strchr(id, ' ')))
				*ch = '\0';
	//	}
	//	else {
	//		if (sscanf(buf, "%s %c %s %d %d", id, &ori, chrom, &b, &l) != 5)
	//			errAbort("error: %s", buf);
	//	}
		lovlp = rovlp = FALSE;
		strcpy(fub, id);
		if (ncbi == FALSE) {
			len = strlen(fub);
//			assert(fub[len - 4] == '#');
			assert(fub[len - 1] == '2');
			fub[len - 1] = '1';
		}
		el = hashLookup(leftReads, fub);
		while (el) {
			len = l;
			pp = (struct pos *)(el->val);
			tr = (struct rbTree *)hashFindVal(homoHash, pp->chr);
			if(tr != NULL){
				lrg = rangeTreeMaxOverlapping(tr, pp->beg, pp->end);
			}
			else{
				lrg = NULL;
			}

			tr = (struct rbTree *)hashFindVal(homoHash, chrom);
			if(tr != NULL){
				rrg = rangeTreeMaxOverlapping(tr, b, b + len - 1);
			}
			else{
				rrg = NULL;
			}
			if (lrg) {
				if (min(lrg->end, pp->end) - max(lrg->start, pp->beg) + 1 > len * D)
					lovlp = TRUE;
			}
			if (rrg) {
					if (min(rrg->end, b + len - 1) - max(rrg->start, b) + 1 > len * D)
					rovlp = TRUE;
			}
			if (lovlp || rovlp) {
				el = hashLookupNext(el);
				continue;
			}
			printf("%s\t%s %d %d %c\t%s %d %d %c\n", 
						fub, pp->chr, pp->beg, pp->end, pp->orient, chrom, b, b+len-1, ori);
			el = hashLookupNext(el);
		}
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 4)
		usage();
	ncbi = optionExists("ncbi");
	optionMustExist("repeatOverlap");
	D = optionDouble("repeatOverlap", 0.5);
	leftReads = homoHash = NULL;
	
	createRangeTreeHash(argv[3]);
	createReadsHash(argv[1]);
	leftRightOvlp(argv[2]);
	
	hashFreeWithVals(&leftReads, freez);
	hashFreeWithVals(&homoHash, rbTreeFree);

	return 0;
}

