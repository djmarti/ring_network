#ifndef _PARSER_H
#define _PARSER_H 1
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "ext_inputs.h"

char *rstrip(char *);
char *lskip(char *);
char *strip_quotes(char *);
char *find_char_or_comment(char *, char);
char *strncpy0(char *dest, const char *, size_t);
int conf_parse(const char *filename,
               int (*handler) (void *usercfg, const char *name,
                               const char *value), void *usercfg);
int extinputs_parse(const char *filename, struct pulse **cfginputs, int N);
int extinputs_parse_2D(const char *filename, struct pulse_2D **cfginputs,
                       int N1, int N2);

#endif
