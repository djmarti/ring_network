#ifndef __CONF_FILE_H_
#define __CONF_FILE_H_ 1
#include <getopt.h>
#include <stdio.h>
#include "parameters.h"
#include "parser.h"

/* conf_file.h -- simple conf file parser

   Based on inih, released under the New BSD license (see LICENSE.txt). Go to the
   project home page for more info:

   http://code.google.com/p/inih/ 



   Parse given conf-style file with name = value pairs
   (whitespace stripped), and comments starting with '#' (number sign). For
   each name=value pair parsed, call handler function with given user pointer
   as well as  name and value (data only valid for duration of handler call).
   Handler should return nonzero on success, zero on error.

   Returns 0 on success, line number of first error on parse error, or -1 on
   file open error.
 */

int process_options(void *, int, char *[]);
/* The next is an abridged version of process_options necessary for
 * updating the config file without parsing all the options */
int process_config_in_options(void *usercfg, int, char *[]);
void usage(int, char *);
#endif
