#include <stdio.h>
#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"

#define AT_RATIO	0.4
#define PL			0.3
#define LEFT		0
#define RIGHT		1
#define AVG			10

void usage() {
	errAbort(
		"singleFusionJunction - summarize junction spanning reads for one candidate fusion.\n"
		"    Usage: singleFusionJunction <knownExonIntron> <output_dir> <output_file> -minSpan=? -minOverlap=?"
	);
}

struct info {
	char *dna;
	char *pos[2], *gene[2];
};

struct segStruct {
	struct segStruct *next;
	char *name;
	boolean poly;
};

static struct optionSpec options[] = {
	{"minSpan", OPTION_INT},
	{"minOverlap", OPTION_INT},
	{NULL, 0},
};

static struct hash *exonHash = NULL, *intronHash = NULL;
static int MINSPAN, MINJ;

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

static void freeInfoSpace(struct info **inf) {
	freeMem((*inf)->dna);
	freeMem((*inf)->pos[LEFT]);
	freeMem((*inf)->pos[RIGHT]);
	freeMem((*inf)->gene[LEFT]);
	freeMem((*inf)->gene[RIGHT]);
	freez(inf);
}

static void createGeneHash(char *geneList) {
	FILE *fp;
	char buf[500], gene[50], chr[50];
	int b, e;
	struct slName *sn;
	struct hashEl *el;
	struct rbTree *tr;

	exonHash = newHash(8);
	intronHash = newHash(8);

	fp = mustOpen(geneList, "r");
	while (fgets(buf, 500, fp)) {
		if (buf[0] == '@') {
			if (sscanf(buf, "@ %s %d %d %s", chr, &b, &e, gene) != 4)
				errAbort("error: %s", buf);
			sn = newSlName(gene);
			el = hashLookup(exonHash, chr);
			if (el == NULL) {
				tr = rangeTreeNew();
				hashAdd(exonHash, chr, tr);
			}
			else
				tr = (struct rbTree *)(el->val);
			rangeTreeAddValHead(tr, b, e, &sn);
		}
		else if (buf[0] == '-') {
			if (sscanf(buf, "- %s %d %d %s", chr, &b, &e, gene) != 4)
				errAbort("error: %s", buf);
			sn = newSlName(gene);
			el = hashLookup(intronHash, chr);
			if (el == NULL) {
				tr = rangeTreeNew();
				hashAdd(intronHash, chr, tr);
			}
			else
				tr = (struct rbTree *)(el->val);
			rangeTreeAddValHead(tr, b, e, &sn);
		}
		else errAbort("error: %s", buf);
	}
	fclose(fp);
}

static boolean slNameNotUnique(struct slName *list, struct slName *fl) {
	struct slName *el;
	boolean afterThis = FALSE, hasSame = FALSE;
	for (el = list; el != NULL; el = el->next) {
		if (el != fl && sameWord(fl->name, el->name)) {
			hasSame = TRUE;
			break;
		}
		if (el == fl)
			afterThis = TRUE;
	}
	if (hasSame == TRUE && afterThis == FALSE)
		return TRUE;
	return FALSE;
}

static boolean maybePolyA(char *dna) {
	int len;
	int i, aa, tt, j;
	char abuf[500], bbuf[500];
	char str[2][500];
	char *ch;
	boolean poly = FALSE;
	strcpy(str[0], dna);
	ch = strchr(dna, ' ');
	++ch;
	strcpy(str[1], ch);
	ch = strchr(str[0], ' ');
	*ch = '\0';
	for (j = LEFT; j <= RIGHT; j++) {
		len = strlen(str[j]);
		for (i = aa = tt = 0; i < len; i++) {
			if (str[j][i] == 'T')
				tt++;
			else if (str[j][i] == 'A')
				aa++;
		}
		for (i = 0; i < len * PL; i++) {
			abuf[i] = 'A';
			bbuf[i] = 'T';
		}
		abuf[i] = bbuf[i] = '\0';
		if ((strstr(str[j], abuf) && aa >= AT_RATIO * len ) || 
				(strstr(str[j], bbuf) && tt >= AT_RATIO * len) )
			poly = TRUE;
	}
	return poly;
}

void findFusions(char *dir, char *out) {
	FILE *fp = NULL, *fd;
	char buf[500], chrom[500], fub[500], ori[2], geneName[500], read[500], cigar[50], dna[500], key[500];
	char str[2][500], chr[2][50];
	char *ch;
	int rlen = 0, i, total, off, k, m, init, last, j, front, gap, end, pos, num, segNum;
	int b[2], e[2], s[2], t[2], len[2], sideSize[2], totalSize[2];
	struct hash *readsHash = NULL, *tmpHash = NULL;
	struct hashEl *el, *gl;
	struct hashCookie cookie;
	struct info *inf;
	struct rbTree *tr;
	struct range *rg;
	struct slName *sn;
	struct segStruct *seg, *head;
	long size;
	struct hash *leftBegEnd, *rightBegEnd;

	sprintf(fub, "%s/reads.fq", dir);
	fd = mustOpen(fub, "r");
	fseek(fd , 0 , SEEK_END);
	size = ftell(fd);
	fclose(fd);
	if (size == 0)
		return;
	
	fd = mustOpen(out, "w");
	sprintf(buf, "%s", dir);
	sprintf(fub, "%s/output/final.sam", buf);
	fp = mustOpen(fub, "r");
	readsHash = hashNew(8);
	tmpHash = hashNew(8);
	j = init = last = 0;
	leftBegEnd = hashNew(4);
	rightBegEnd = hashNew(4);
	while (fgets(buf, 500, fp)) {
		if (buf[0] == '@')
			continue;
		if (sscanf(buf, "%s %*d %s %d %*s %s %*s %*s %*s %s %*s", read, chrom, &pos, cigar, dna) != 5) 
			errAbort("error: %s", buf);
		if (rlen == 0) {
			rlen = strlen(read);
			rlen += 6;
		}
		if (strlen(cigar) > 3) {
			if (sscanf(cigar, "%dM%dN%dM", &front, &gap, &end) != 3)
				errAbort("cigar: %s", cigar);
			if (front < MINJ || end < MINJ)
				continue;
			if (sscanf(chrom, "%[^:]:%d-%d+%[^:]:%d-%d[%c%c](%d)", 
												chr[LEFT], &(b[LEFT]), &(e[LEFT]), 
												chr[RIGHT], &(b[RIGHT]), &(e[RIGHT]), &(ori[LEFT]), &(ori[RIGHT]), &num) != 9)
				errAbort("chr: %s", chrom);
			len[LEFT] = e[LEFT] - b[LEFT] + 1;
			len[RIGHT] = e[RIGHT] - b[RIGHT] + 1;
			if (pos < len[LEFT] && pos + front + gap > len[LEFT]) {
				
				AllocVar(inf);
				if (ori[LEFT] == '+')
					sprintf(fub, "%s:%d-%d", chr[LEFT], b[LEFT] + pos - 1, b[LEFT] + pos + front - 2);
				else
					sprintf(fub, "%s:%d-%d", chr[LEFT], e[LEFT] - (pos + front - 2), e[LEFT] - (pos - 1));
				inf->pos[LEFT] = cloneString(fub);
				off = pos + front + gap - len[LEFT] - 1;
				if (ori[RIGHT] == '+')
					sprintf(fub, "%s:%d-%d", chr[RIGHT], b[RIGHT] + off, b[RIGHT] + off + end - 1);
				else
					sprintf(fub, "%s:%d-%d", chr[RIGHT], e[RIGHT] - off - end + 1, e[RIGHT] - off );//+1
				inf->pos[RIGHT] = cloneString(fub);
				
				for (i = LEFT; i <= RIGHT; i++) {
					sscanf(inf->pos[i], "%*[^:]:%d-%d", &s[i], &t[i]);
					gl = hashLookup(exonHash, chr[i]);
					tr = (struct rbTree *)(gl->val);
					rg = rangeTreeMaxOverlapping(tr, s[i], t[i]);
					if (rg) {
						for (sn = (struct slName *)(rg->val); sn; sn = sn->next) {
							if (!slNameNotUnique((struct slName *)(rg->val), sn))
								sprintf(geneName, " %s", sn->name);
						}
						inf->gene[i] = cloneString(geneName);
					}
					else {
					/*
						gl = hashLookup(intronHash, chr[i]);
						tr = (struct rbTree *)(gl->val);
						rg = rangeTreeMaxOverlapping(tr, s[i], t[i]);
						if (rg) {
							for (sn = (struct slName *)(rg->val); sn; sn = sn->next) {
								if (!slNameNotUnique((struct slName *)(rg->val), sn))
									sprintf(geneName, " %s", sn->name);
								inf->gene[i] = cloneString(geneName);
							}
						}
						else inf->gene[i] = NULL;
					*/
						inf->gene[i] = NULL;
					}
				}

				if (inf->gene[LEFT] == NULL || inf->gene[RIGHT] == NULL) {
					freeInfoSpace(&inf);
					continue;
				}
				for (k = 0; k < front; k++)
					fub[k] = dna[k];
				fub[k++] = ' ';
				m = strlen(dna);
				for (; k < m+1; k++) 
					fub[k] = dna[k-1];
				fub[k] = '\0';
				inf->dna = cloneString(fub);

/*				if (hashLookup(leftBegEnd, inf->pos[LEFT]))
					continue;

				else
					hashStoreName(leftBegEnd, inf->pos[LEFT]);		
				if (hashLookup(rightBegEnd, inf->pos[RIGHT]))
					continue;
	
				else
					hashStoreName(rightBegEnd, inf->pos[RIGHT]);*/
				if (ori[LEFT] == '+' && ori[RIGHT] == '+')
					sprintf(fub, "%d %d", t[LEFT], s[RIGHT]);
				else if (ori[LEFT] == '+' && ori[RIGHT] == '-')
					sprintf(fub, "%d %d", t[LEFT], t[RIGHT]);
				else if (ori[LEFT] == '-' && ori[RIGHT] == '+')
					sprintf(fub, "%d %d", s[LEFT], s[RIGHT]);
				else if (ori[LEFT] == '-' && ori[RIGHT] == '-')
					sprintf(fub, "%d %d", s[LEFT], t[RIGHT]);
		
				sprintf(key, "%s %s", read, fub);
				if (hashLookup(readsHash, key)) {
					freeInfoSpace(&inf);
					continue;
				}
				hashAdd(readsHash, key, inf);
				
				AllocVar(seg);
				seg->name = cloneString(read);
				seg->poly = maybePolyA(inf->dna);
				el = hashLookup(tmpHash, fub);
				if (el) {
					head = (struct segStruct *)(el->val);
					slAddHead(&head, seg);
					el->val = head;
				}
				else {
					hashAdd(tmpHash, fub, seg);
					
				}
				++j;
			}
		}
	}
	fclose(fp);
	hashFree(&leftBegEnd);
	hashFree(&rightBegEnd);
	if (j > 0) {
		total = 0;
		cookie = hashFirst(tmpHash);
		while((el = hashNext(&cookie))) {
			i = 0;
			init = last = totalSize[0] = totalSize[1] = sideSize[0] = sideSize[1] = 0;
			for (seg = (struct segStruct *)(el->val); seg; seg = seg->next) {
				if (seg->poly)
					++i;
				sprintf(key, "%s %s", seg->name, el->name);
				inf = (struct info *)hashMustFindVal(readsHash, key);
				strcpy(str[LEFT], inf->dna);
				ch = strchr(str[LEFT], ' ');
				*ch = '\0';
				++ch;
				strcpy(str[RIGHT], ch);
				if (strlen(str[LEFT]) <= (unsigned)MINJ)
					++sideSize[LEFT];
				if (strlen(str[RIGHT]) <= (unsigned)MINJ)
					++sideSize[RIGHT];
				totalSize[LEFT] += strlen(str[LEFT]);
				totalSize[RIGHT] += strlen(str[RIGHT]);
				front = ch - str[LEFT] - 1;
				end = strlen(inf->dna) - front - 1;
				init = max(init, front);
				last = max(last, end);
			}
			segNum = slCount((struct segStruct *)(el->val));

			if (i >= 0.5 * segNum || segNum < MINSPAN || 
					sideSize[LEFT] == segNum || sideSize[RIGHT] == segNum)
				continue;

			if (totalSize[LEFT]/segNum <= AVG || totalSize[RIGHT]/segNum <= AVG)
				continue;
			if (segNum == 1){
				seg = (struct segStruct *)(el->val);
				sprintf(key, "%s %s", seg->name, el->name);
                                inf = (struct info *)hashMustFindVal(readsHash, key);
                                strcpy(str[LEFT], inf->dna);
                                ch = strchr(str[LEFT], ' ');
                                *ch = '\0';
                                ++ch;
                                strcpy(str[RIGHT], ch);
				if(strlen(str[LEFT]) < 15 || strlen(str[RIGHT]) < 15)
					continue;
			}
			if (total != 0)
				fprintf(fd, "----------------\n");
			for (seg = (struct segStruct *)(el->val); seg; seg = seg->next) {
				
				sprintf(key, "%s %s", seg->name, el->name);
				inf = (struct info *)hashMustFindVal(readsHash, key);
				if(! strcmp(inf->gene[LEFT],inf->gene[RIGHT]))
					continue;
				++total;
				if (total == 1) 
					fprintf(fd, "# Fusion: %s %s\n", dir, chrom);
				
//				inf = (struct info *)hashMustFindVal(readsHash, key);
				strcpy(str[LEFT], inf->dna);
				ch = strchr(str[LEFT], ' ');
				*ch = '\0';
				++ch;
				strcpy(str[RIGHT], ch);
				fprintf(fd, "-> %-*s%*s %-*s  %s %s ", rlen, seg->name, init, str[LEFT], last, str[RIGHT], inf->pos[LEFT], inf->pos[RIGHT]);
				fprintf(fd, "%s x%s\n", inf->gene[LEFT], inf->gene[RIGHT]);
			}
		}
		if (total > 0) 
			fprintf(fd, "# Total # of Junction Spanning Reads: %d\n\n", total);
	}
	fclose(fd);
	hashFreeWithVals(&tmpHash, slFreeList);
	hashFreeWithVals(&readsHash, freeInfoSpace);
	hashFreeWithVals(&exonHash, rbTreeFree);
	hashFreeWithVals(&intronHash, rbTreeFree);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 4)
		usage();
	optionMustExist("minOverlap");
	MINJ = optionInt("minOverlap", 6);
	optionMustExist("minSpan");
	MINSPAN = optionInt("minSpan", 2);
	createGeneHash(argv[1]);
	findFusions(argv[2], argv[3]);
	return 0;
}

