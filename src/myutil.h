#ifndef _UTIL_
#define _UTIL_

#define NIL -1
#define HSIZE	12
#define BIG 9999999

void swapInt(int *a, int *b);

void do_cmd(const char *fmt, ...);

int overlap(int b1, int e1, int b2, int e2);

#endif
