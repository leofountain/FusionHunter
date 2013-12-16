#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"
#include "myutil.h"
#include <omp.h>


void usage() {
	errAbort(
		"allFusionJunction - summarize all the junctions.\n"
		"    Usage: allFusionJunction <region_file> <knownAnnotation> <output_dir> <output-file> -minSpan=? -minOverlap=?"
	);
}

static struct optionSpec options[] = {
	{"minSpan", OPTION_INT},
	{"minOverlap", OPTION_INT},
	{NULL, 0},
};

static int MINSPAN, MINJ;

void findFusions(char *regions, char *geneList, char *dir, char *fileName, char *bin) {
	FILE *fp = NULL;
	char buf[500], dirName[500];
	int i, total;

	fp = mustOpen(regions, "r");
	total = 0;
	while(fgets(buf, 500, fp))
		++total;
	fclose(fp);

	#pragma omp parallel for private(i, dirName)
	for (i = 1; i <= total; i++) {
		sprintf(dirName, "%s/R%d", dir, i);
		do_cmd("%s/singleFusionJunction -minOverlap=%d -minSpan=%d %s %s %s/output/out.junctions",bin, MINJ, MINSPAN, geneList, dirName, dirName);
		do_cmd("%s/ifJunctionAnnotated.pl %s/output/final.sam %s/output/out.junctions > %s/output/out.fusions",bin, dirName,dirName,dirName);
		if(strcmp(regions,"RG"))do_cmd("%s/fusionInterpret.pl %s/output/out.fusions Z_known.exon.intron %s/output",bin,dirName,dirName);
	}

	for (i = 1; i <= total; i++) {
		sprintf(dirName, "%s/R%d", dir, i);
		if (i == 1)
			do_cmd("cat %s/output/out.fusions > %s", dirName, fileName);
		else
			do_cmd("cat %s/output/out.fusions >> %s", dirName, fileName);
	}
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 6)
		usage();
	optionMustExist("minOverlap");
	MINJ = optionInt("minOverlap", 6);
	optionMustExist("minSpan");
	MINSPAN = optionInt("minSpan", 2);
	findFusions(argv[1], argv[2], argv[3], argv[4], argv[5]);
	return 0;
}

