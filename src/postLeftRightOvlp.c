#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"
#include "myutil.h"

struct aliStruct {
	struct aliStruct *next;
	char chr1[10], chr2[10];
	int beg1, beg2, end1, end2, id;
};

void usage() {
	errAbort(
		"postLeftRightOvlp - retain paired-end reads whose left and right do not share homology\n"
		"    Usage: postLeftRightOvlp <reads> <homology> -chainDist=? -chainOverlap=?"
	);
}

static struct optionSpec options[] = {
	{"chainDist", OPTION_INT},
	{"chainOverlap", OPTION_DOUBLE},
	{NULL, 0},
};

struct hash *aliHash = NULL;
static int FEW;
static double D;

static struct range *rangeTreeAddValHead(struct rbTree *tree, int start, int end, struct aliStruct **newVal) {
	struct range *r, *existing;
	struct aliStruct *head;
	AllocVar(r);
	r->start = start;
	r->end = end;
	r->val = *newVal;
	while ((existing = rbTreeRemove(tree, r))) {
    r->start = min(r->start, existing->start);
    r->end = max(r->end, existing->end);
		head = (struct aliStruct *)(existing->val);
    slAddHead(&head, *newVal);
		r->val = head;
	}
	rbTreeAdd(tree, r);
	return r;
}

static void readSelfAlignment(char *fileName) {
	FILE *fp;
	char buf[500], chrom1[10], chrom2[10];
	int b1, b2, e1, e2, id;
	struct aliStruct *ali;
	struct hashEl *el;
	struct rbTree *tr;

	aliHash = newHash(8);
	fp = mustOpen(fileName, "r");
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%s %d %d %s %d %d %*c [%d]", chrom1, &b1, &e1, chrom2, &b2, &e2, &id) != 7)
			errAbort("error: %s", buf);
		AllocVar(ali);
		strcpy(ali->chr1, chrom1);
		strcpy(ali->chr2, chrom2);
		ali->beg1 = b1;
		ali->end1 = e1;
		ali->beg2 = b2;
		ali->end2 = e2;
		ali->id = id;
		el = hashLookup(aliHash, chrom1);
		if (el == NULL) {
			tr = rangeTreeNew();
			hashAdd(aliHash, chrom1, tr);
		}
		else 
			tr = (struct rbTree *)(el->val);
		rangeTreeAddValHead(tr, b1, e1, &ali);
	}
	fclose(fp);
}

/*
static void printHash() {
	struct hashEl *el;
	struct hashCookie cookie;
	struct rbTree *tr;
	struct range *rg;
	struct aliStruct *ali;
	
	cookie = hashFirst(aliHash);
	while ((el = hashNext(&cookie))) {
		tr = (struct rbTree *)(el->val);
		for (rg = rangeTreeList(tr); rg; rg = rg->next) {
			fprintf(stderr, "# %s %d %d\n", el->name, rg->start, rg->end);
			for (ali = (struct aliStruct *)(rg->val); ali; ali = ali->next) {
				fprintf(stderr, "-> %s %d %d %s %d %d [%d]\n", ali->chr1, ali->beg1, ali->end1, ali->chr2, ali->beg2, ali->end2, ali->id);
			}
		}
	}
} */

void postLeftRightOvlp(char *readsFile) {
	FILE *fp;
	char buf[500], chrom1[10], chrom2[10];
	int b1, b2, e1, e2, len1, len2;
	struct aliStruct *ali;
	struct hashEl *el;
	struct rbTree *tr;
	struct range *rg;
	boolean left, right;

	fp = mustOpen(readsFile, "r");
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%*s %s %d %d %*c %s %d %d %*c", chrom1, &b1, &e1, chrom2, &b2, &e2) != 6)
			errAbort("error: %s", buf);
		if (sameString(chrom1, "chrM") || sameString(chrom1, "chrY") 
				|| sameString(chrom2, "chrM") || sameString(chrom2, "chrY"))
			continue;
		left = right = FALSE;
		el = hashLookup(aliHash, chrom1);
		if(el == NULL){
			printf("%s", buf);
			continue;
		}

		//assert(el);
		tr = (struct rbTree *)(el->val);
		for (rg = rangeTreeAllOverlapping(tr, b1, e1); rg; rg = rg->next) {
			for (ali = (struct aliStruct *)(rg->val); ali; ali = ali->next) {
				if (sameString(ali->chr2, chrom2)) {
					len1 = e1 - b1 + 1;
					len2 = e2 - b2 + 1;
					if (overlap(b1, e1, max(0, ali->beg1 - FEW), ali->end1 + FEW) > len1 * D && 
							overlap(b2, e2, max(0, ali->beg2 - FEW), ali->end2 + FEW) > len2 * D)  {
						break;
					}
				}
			}
			if (ali) {
				left = TRUE;
				break;
			}
		}
		if (left == FALSE) {
			el = hashLookup(aliHash, chrom2);
			//assert(el);
			if(el == NULL){
				printf("%s", buf);
				continue;
			}

			tr = (struct rbTree *)(el->val);
			for (rg = rangeTreeAllOverlapping(tr, b2, e2); rg; rg = rg->next) {
				for (ali = (struct aliStruct *)(rg->val); ali; ali = ali->next) {
					if (sameString(ali->chr2, chrom1)) {
						len1 = e1 - b1 + 1;
						len2 = e2 - b2 + 1;
						if (overlap(b1, e1, max(0, ali->beg2 - FEW), ali->end2 + FEW) > len1 * D && 
								overlap(b2, e2, max(0, ali->beg1 - FEW), ali->end1 + FEW) > len2 * D)  {
							break;
						}
					}
				}
				if (ali) {
					right = TRUE;
					break;
				}
			}
		}
		if (left == FALSE && right == FALSE) 
			printf("%s", buf);
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 3)
		usage();
	optionMustExist("chainOverlap");
	optionMustExist("chainDist");
	FEW = optionInt("chainDist", 100000);
	D = optionDouble("chainOverlap", 0.5);
	readSelfAlignment(argv[2]);
	postLeftRightOvlp(argv[1]);
	
	hashFreeWithVals(&aliHash, rbTreeFree);
	
	return 0;
}

