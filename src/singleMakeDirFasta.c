#include "common.h"
#include "errabort.h"
#include "options.h"
#include "fa.h"
#include "dnaseq.h"
#include "hash.h"
#include "myutil.h"

void usage() {
	errAbort(
		"singleMakeDirFasta - make directories for one the regions and create fasta files\n"
		"    Usage: makeDirFasta <regions> <hg18_fasta_file> <dir> <num>"
	);
}

static struct optionSpec options[] = {
	{NULL, 0},
};

static char *cloneStringZExt(const char *s, int size, int copySize)
/* Make a zero terminated copy of string in memory */
{
	char *d = needMem(copySize+1);
	copySize = min(size,copySize);
	memcpy(d, s, copySize);
	d[copySize] = 0; 
	return d; 
}

static void makeDirFasta(char *regionsFile, char *hg18FastaFile, char *dir, int num) {
	FILE *fp, *sq;
	char buf[500], dirName[500], seqName[500], chr1[500], chr2[500];
	int b1, e1, b2, e2, i, len;
	char ori1, ori2;
	struct hash *seqHash = NULL;
	struct dnaSeq *seq1, *seq2;
	struct stat st;
	DNA *s1, *s2;

	seqHash = faReadAllIntoHash(hg18FastaFile, dnaUpper);
	if (stat(dir, &st) != 0)
		do_cmd("mkdir %s", dir);

	fp = mustOpen(regionsFile, "r");
	i = 0;
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%[^:]:%d-%d %[^:]:%d-%d [%c %c]", chr1, &b1, &e1, chr2, &b2, &e2, &ori1, &ori2) != 8)
			errAbort("error: %s", buf);
		++i;
		if (i != num) 
			continue;
		sprintf(dirName, "%s/R%d", dir, i);
		if (stat(dirName, &st) != 0)
			do_cmd("mkdir %s", dir);
		sprintf(seqName, "%s/ref.fa", dirName);
		sq = mustOpen(seqName, "w");
		fprintf(sq, ">%s:%d-%d+%s:%d-%d[%c%c]\n", chr1, b1, e1, chr2, b2, e2, ori1, ori2);
		seq1 = (struct dnaSeq *)hashFindVal(seqHash, chr1);
		assert(e1 <= seq1->size);
		len = e1 - b1 + 1;
		if (ori1 == '-') {
			s1 = cloneStringZExt(seq1->dna + b1 - 1, len, len+1);
			reverseComplement(s1, len);
			writeSeqWithBreaks(sq, s1, len, 80);
			freeMem(s1);
		}
		else
			writeSeqWithBreaks(sq, seq1->dna + b1 - 1, e1 - b1 + 1, 80);
		seq2 = (struct dnaSeq *)hashFindVal(seqHash, chr2);
		assert(e2 <= seq2->size);
		len = e2 - b2 + 1;
		if (ori2 == '-') {
			s2 = cloneStringZExt(seq2->dna + b2 - 1, len, len+1);
			reverseComplement(s2, len);
			writeSeqWithBreaks(sq, s2, len, 80);
			freeMem(s2);
		}
		else
			writeSeqWithBreaks(sq, seq2->dna + b2 - 1, e2 - b2 + 1, 80);
		fclose(sq);
	}
	fclose(fp);
	//FIXME: free space
} 

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 5)
		usage();
	makeDirFasta(argv[1], argv[2], argv[3], atoi(argv[4]));
	return 0;
}

