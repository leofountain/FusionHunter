#include "common.h"
#include "errabort.h"
#include "options.h"
#include "hash.h"
#include "rangeTree.h"
#include "obscure.h"



void usage() {
	errAbort(
		"regionPairsList - print list of region association before running the rest\n"
		"    Usage: regionPairsList <regions> -minPair=?"
	);
}

static struct optionSpec options[] = {
	{"minPair", OPTION_INT},
	{NULL, 0},
};

static int MINPAIR;

// based on illumina paired-end reads
// FR						AB 
// -+	=>	0x00	-- 
// --	=>	0x01	-+
// ++	=>	0x10	+-
// +-	=>	0x11	++
static void mapOrient(int i, char *ori) {
	if (i == 0) sprintf(ori, "[- -]");
	if (i == 1) sprintf(ori, "[- +]");
	if (i == 2) sprintf(ori, "[+ -]");
	if (i == 3) sprintf(ori, "[+ +]");
}

static void regionPairsList(char *regionsFile) {
	FILE *fp;
	char buf[500], str1[500], str2[500], orient[10];
	int i;
	int od[4];

	fp = mustOpen(regionsFile, "r");
	while (fgets(buf, 500, fp)) {
		if (sscanf(buf, "%[^ ] %s %d %d %d %d", str1, str2, &(od[0]), &(od[1]), &(od[2]), &(od[3])) != 6)
			errAbort("error: %s", buf);
		for (i = 0; i < 4; i++) {
			if (od[i] >= MINPAIR) {
				mapOrient(i, orient);
				if (sameString(orient, "[- -]"))
					printf("%s\t%s\t[+ +]\t(%d)\n", str2, str1, od[i]);
				else
					printf("%s\t%s\t%s\t(%d)\n", str1, str2, orient, od[i]);
			}
		}
	}
	fclose(fp);
} 

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	optionMustExist("minPair");
	MINPAIR = optionInt("minPair", 4);
	regionPairsList(argv[1]);
	return 0;
}

