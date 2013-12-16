#include "common.h"
#include "errabort.h"
#include "myutil.h"

void do_cmd(const char *fmt, ...) {
	char buf[10000];
	va_list ap;
	va_start(ap, fmt);
	(void)vsprintf(buf, fmt, ap);
	if (system(buf) != 0)
		errAbort("Command '%s' failed", buf);
	va_end(ap);
}

void swapInt(int *a, int *b) {
	int t;
	t = *a;
	*a = *b;
	*b = t;
}

int overlap(int b1, int e1, int b2, int e2) {
	int op = 0;
	if (b1 <= e2 && e1 >= b2)
		op = min(e1, e2) - max(b1, b2) + 1;
	assert(op >= 0);
	return op;
}

