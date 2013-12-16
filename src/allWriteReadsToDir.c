#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

#define D 0.2
#define FEW 10000

void usage() {
	errAbort(
		"allWriteReadsToDir - get reads for each region\n"
		"    Usage: allWriteReadsToDir <regions> <original-left> <original-right> <half-left> <half-right> <full-left> <full-right> <dir>"
	);
}

static struct optionSpec options[] = {
	{"ncbi", OPTION_BOOLEAN},
	{NULL, 0},
};

struct hash *aliHash = NULL, *readsHash = NULL;
static boolean ncbi = FALSE;

static struct range *rangeTreeAddValHead(struct rbTree *tree, int start, int end, struct slName **newVal) {
	struct range *r, *existing;
	struct slName *head;
	AllocVar(r);
	r->start = start;
	r->end = end;
	r->val = *newVal;
	while ((existing = rbTreeRemove(tree, r))) {
    		r->start = min(r->start, existing->start);
		r->end = max(r->end, existing->end);
		head = (struct slName *)(existing->val);
   		slAddHead(&head, *newVal);
		r->val = head;
	}
	rbTreeAdd(tree, r);
	return r;
}

static void readMappingResults(char **argv) {
	FILE *fp;
	char buf[500], chr[50], id[500];
	int beg, i, j, len;
	char *ch;
	struct slName *ali;
	struct hashEl *el;
	struct rbTree *tr;

	aliHash = newHash(8);
	
	for (i = 4; i <= 5; i++) {
		fp = mustOpen(argv[i], "r");
		while (fgets(buf, 500, fp)) {
			
			if (ncbi && i >= 6 && i <= 7) {
				if (sscanf(buf, "%[^\t]\t%*c %s %d %d", id, chr, &beg, &len) != 4)
					errAbort("error: %s", buf);
				if ((ch = strchr(id, ' ')))
					*ch = '\0';
				if (i >= 6 && i <= 7) {
					j = strlen(id);
					sprintf(id+j, "/%d", i-5);
				}
			}
			
			else {	
				if (sscanf(buf, "%[^\t]\t%*c %s %d %d", id, chr, &beg, &len) != 4){
					errAbort("error: %s", buf);
				}
					if ((ch = strchr(id, ' ')))
                                        	*ch = '\0';
				
			}
			ali = newSlName(id);
			el = hashLookup(aliHash, chr);
			if (el == NULL) {
				tr = rangeTreeNew();
				hashAdd(aliHash, chr, tr);
			}
			else 
				tr = (struct rbTree *)(el->val);
			rangeTreeAddValHead(tr, beg, beg + len - 1, &ali);
		}
		fclose(fp);
	}
}
void getRange(char *regionFile){
                FILE *fp;
        char buf[500],chr[50];
        char str[2][500];
        int i, b, e;
        struct hashEl *el;
        struct rbTree *tr;
        fp = mustOpen(regionFile, "r");
        aliHash = newHash(8);
        while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%[^\t]\t%[^\t]\t%*s", str[0], str[1]) != 2)
                        errAbort("error: %s", buf);
		for (i = 0; i < 2; i++) {
                        if (sscanf(str[i], "%[^:]:%d-%d", chr, &b, &e) != 3)
                                errAbort("error: %s", str[i]);
     //   aliHash = newHash(8);
                	el = hashLookup(aliHash, chr);
                	if (el == NULL) {
                                tr = rangeTreeNew();
                                hashAdd(aliHash, chr, tr);
                        }
                	else
                        	tr = (struct rbTree *)(el->val);
	                rangeTreeAdd(tr, b, e);
        	}
	}
//      printf("range\n");
}

static void	createReadsHash(char **argv) {
	int i;
	FILE *fp;
	char buf[500], fub[500], dna[500], qua[500], id[50];
	char *str;
	readsHash = newHash(16);
	for (i = 2; i <= 3; i++) {
		fp = mustOpen(argv[i], "r");
		for (;;) {
			if (fgets(buf, 500, fp)) {
				if(strlen(buf) == 0){break;}
				if (ncbi) {
					sscanf(buf, "@%s %*s", id);
					sprintf(fub, "%s/%d", id, i - 1);
				}
				else {
					sscanf(buf, "@%s %*s", id);
					strcpy(fub, id);
				}
				if (!fgets(buf, 500, fp)){
					break;
					errAbort("error: %s", argv[i]);
	}
				sscanf(buf, "%s", dna);
				if (!fgets(buf, 500, fp))
					errAbort("error: %s", argv[i]);
				if (!fgets(buf, 500, fp))
					errAbort("error: %s", argv[i]);
				sscanf(buf, "%s", qua);
				sprintf(buf, "%s %s", dna, qua);
				str = cloneString(buf);
				hashAdd(readsHash, fub, str);
			}
			else break;
		}
		fclose(fp);
	}
}


void allWriteReadsToDir(char *regionFile, char *dir) {
	FILE *fp, *rd;
	char buf[500], readName[500], fileName[500], chr[50], fub[500];
	char str[2][500];
	char *readStr, *ch;
	int i, b, e, j, k;
	struct slName *ali;
	struct hashEl *el;
	struct rbTree *tr;
	struct range *rg;
	struct hash *localHash = NULL;

	fp = mustOpen(regionFile, "r");
	j = 0;
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%[^\t]\t%[^\t]\t%*s", str[0], str[1]) != 2)
			errAbort("error: %s", buf);
		++j;
		sprintf(fileName, "%s/R%d/reads.fq", dir, j);
		rd = mustOpen(fileName, "w");
		localHash = hashNew(8);
		for (i = 0; i < 2; i++) {
			if (sscanf(str[i], "%[^:]:%d-%d", chr, &b, &e) != 3)
				errAbort("error: %s", str[i]);
			el = hashLookup(aliHash, chr);
			tr = (struct rbTree *)(el->val);
			for (rg = rangeTreeAllOverlapping(tr, b, e); rg; rg = rg->next) {
				for (ali = (struct slName *)(rg->val); ali; ali = ali->next) {
					if (hashLookup(localHash, ali->name))
						continue;
					hashStoreName(localHash, ali->name);
					readStr = (char *)hashFindVal(readsHash, ali->name);
					if(readStr == NULL)
						continue;
					//assert(readStr);
					strcpy(fub, readStr);
					ch = strchr(fub, ' ');
					*ch = '\0';
					fprintf(rd, "@%s\n", ali->name);
					fprintf(rd, "%s\n", fub);
					++ch;
					fprintf(rd, "+%s\n", ali->name);
					fprintf(rd, "%s\n", ch);
         				strcpy(readName, ali->name);
					k = strlen(readName);
					/*
					if (readName[k-1] == '1')
						readName[k-1] = '2';
					else if (readName[k-1] == '2')
						readName[k-1] = '1';
					else
						errAbort("read identifier error: %s", readName);
						
					if (hashLookup(localHash, readName))
						continue;
					hashStoreName(localHash, readName);
					readStr = (char *)hashFindVal(readsHash, readName);

					assert(readStr);
					strcpy(fub, readStr);
					ch = strchr(fub, ' ');
					*ch = '\0';
					fprintf(rd, "@%s\n", readName);
					fprintf(rd, "%s\n", fub);
					++ch;
					fprintf(rd, "+%s\n", readName);
					fprintf(rd, "%s\n", ch); 
					*/
				}
			}
		}
		hashFree(&localHash);
		fclose(rd);
	}
	fclose(fp);
	hashFreeWithVals(&readsHash, freez);
	hashFreeWithVals(&aliHash, rbTreeFree);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	//if (argc != 9)
		//usage();
	ncbi = optionExists("ncbi");
	getRange(argv[1]);
	readMappingResults(argv);
	createReadsHash(argv);
	allWriteReadsToDir(argv[1], argv[6]);
	return 0;
}

