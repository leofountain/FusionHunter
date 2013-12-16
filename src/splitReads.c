#include "common.h"
#include "errabort.h"
#include "options.h"

void usage() {
	errAbort(
		"splitReads solexa-reads -r=? -s=? [-ncbi]\n"
		"   r is the length of the reads' part, e.g. 25\n"
		"   s is the side (0:left, 1:right)"
	);
}

static struct optionSpec options[] = {
	{"r", OPTION_INT},
	{"ncbi", OPTION_BOOLEAN},
	{"s", OPTION_INT},
	{NULL, 0},
};

static int part = 0, side = 0;
static boolean ncbi = FALSE;

void splitReads(char *originalFile) {
	FILE *fp;
	char buf[5000], head1[5000], head2[5000];
	char dna1[5000], dna2[5000], qua1[5000], qua2[5000];
	int len, low, high, r;
	char *ch;

	fp = mustOpen(originalFile, "r");
	r = part;

	for (;;) {
		if (fgets(buf, 5000, fp)) {

			if(buf[0] != '@')
				break;
			sscanf(buf, "%[^\n]", head1);
			
			if (!fgets(buf, 5000, fp))
				errAbort("%s: file error.", originalFile);
			sscanf(buf, "%s", dna1);
			
			if (!fgets(buf, 5000, fp))
				errAbort("%s: file error.", originalFile);
			assert(buf[0] == '+');
			sscanf(buf, "%[^\n]", head2);
			
			if (!fgets(buf, 5000, fp))
				errAbort("%s: file error.", originalFile);
			sscanf(buf, "%s", qua1);
		}
		else break;
		
		low = r;
		len = strlen(dna1);
		if (part > 0.75*len)
			errAbort("r too long (r=%d, readLength=%d)", part, len);
		high = len - r;	//XXX change later to consider splitting a read into more than two parts
		strcpy(dna2, dna1 + high);
		strcpy(qua2, qua1 + high);
		dna1[low] = qua1[low] = '\0';
		
		if (ncbi) {
			if ((ch = strchr(head1, ' '))) {
				*ch ='/';
				sprintf(ch+1, "%d", side+1);
			}
			else {
				ch = head1 + strlen(head1);
				*ch = '/';
				sprintf(ch+1, "%d", side+1);
			}
		}
		printf("%s\n", head1);
		printf("%s\n", dna1);
		
		if (ncbi) {
			if ((ch = strchr(head2, ' '))) {
				*ch = '/';
				sprintf(ch+1, "%d", side+1);
			}
			else {
				ch = head2 + strlen(head2);
				*ch = '/';
				sprintf(ch+1, "%d", side+1);
			}
		}
		printf("%s\n", head2);
		printf("%s\n", qua1);
			
		printf("%s\n", head1);
		printf("%s\n", dna2);
		
		printf("%s\n", head2);
		printf("%s\n", qua2);
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	optionInit(&argc, argv, options);
	if (argc != 2)
		usage();
	optionMustExist("r");
	part = optionInt("r", 25);
	optionMustExist("s");
	side = optionInt("s", 0);
	if (side != 0 && side != 1)
		errAbort("s must be 0 or 1");
	ncbi = optionExists("ncbi");
	splitReads(argv[1]);
	return 0;
}

