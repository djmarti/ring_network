#include "parser.h"

const int MAX_LINE = 200;

/* Auxiliary functions */

/* Strip whitespace chars off end of given string, in place. Return s.  */
char *rstrip(char *s)
{
    char *p = s + strlen(s);
    while (p > s && isspace(*--p))
        *p = '\0';
    return s;
}

/* Return pointer to first non-whitespace char in given string.  */
char *lskip(char *s)
{
    while (*s && isspace(*s))
        s++;
    return (char *) s;
}

/* Return pointer to the unquoted string */
char *strip_quotes(char *s)
{
    char *p = s + strlen(s);
    if (*--p == '"' && *s == '"') {
        s++;
        *p = '\0';
    }
    return s;
}

/* Return pointer to first char c or '#' in given string, or pointer to 
 * null at end of string if neither found.  */
char *find_char_or_comment(char *s, char c)
{
    while (*s && *s != c && *s != '#')
        s++;
    return (char *) s;
}

/* Version of strncpy that ensures dest (size bytes) is null-terminated. */
char *strncpy0(char *dest, const char *src, size_t size)
{
    strncpy(dest, src, size);
    dest[size - 1] = '\0';
    return dest;
}


int conf_parse(const char *filename,
               int (*handler) (void *, const char *, const char *),
               void *usercfg)
{
    FILE *fin;
    char buf[MAX_LINE];

    char *start;
    char *end;
    char *name;
    char *value;

    int error = 0;

    fin = fopen(filename, "r");
    if (!fin)
        return -1;

    int lineno = 0;

    /* Scan through file line by line */
    while (fgets(buf, sizeof(buf), fin) != NULL) {
        lineno++;
        start = lskip(rstrip(buf));     /* chop whites */

        if (*start && *start != '#') {
            /* Not a comment, must be a name = value pair  */
            end = find_char_or_comment(start, '=');
            if (*end == '=') {
                *end = '\0';
                name = rstrip(start);
                value = lskip(end + 1);
                end = find_char_or_comment(value, '#');
                if (*end == '#')
                    *end = '\0';
                value = rstrip(value);
                value = strip_quotes(value);
                /* Valid name=value pair found, call handler  */
                if (handler(usercfg, name, value) && !error)
                    error = lineno;
            } else if (!error) {
                /* No '=' found on name=value line  */
                error = lineno;
            }
        }
    }
    fclose(fin);
    return error;
}

int extinputs_parse(const char *filename, struct pulse **cfginputs, int N)
{
    FILE *fin;
    char buf[MAX_LINE];

    char *start;
    char *end;
    char *values;

    int error = 0;

    fin = fopen(filename, "r");
    if (!fin)
        return -1;

    int lineno = 0;

    double t, q, i, d, m;       /* the values to be read */

    /* Scan through file line by line */
    while (fgets(buf, sizeof(buf), fin) != NULL) {
        lineno++;
        start = lskip(rstrip(buf));     /* chop whites */

        if (*start && *start != '#') {
            /* Remove possible comments at the end */
            end = find_char_or_comment(start, '#');
            if (*end == '#')
                *end = '\0';
            values = rstrip(start);
            /* Not a comment, must be a line containing 4 numbers */
            /* Scan and check if the inputs is correctly formatted */
            if (sscanf(values, "%lf %lf %lf %lf %lf", &t, &d, &q, &i, &m) == 5) {
                if (t >= 0 && d > 0 && m >= 0) {
                    *cfginputs = append_pulse(*cfginputs, new_pulse(N, t, d, q, i, m));
                } else {
                    printf("Durations and time onsets should be positive.\n");
                    printf("Check configuration file for external inputs.\n");
                    error = lineno;
                }
            } else if (!error) {
                /* No '=' found on name=value line  */
                error = lineno;
            }
        }
    }
    fclose(fin);
    return error;
}

int extinputs_parse_2D(const char *filename, struct pulse_2D **cfginputs,
                       int N1, int N2)
{
    FILE *fin;
    char buf[MAX_LINE];

    char *start;
    char *end;
    char *values;

    int error = 0;

    fin = fopen(filename, "r");
    if (!fin)
        return -1;

    int lineno = 0;

    double t, q1, q2, i, d, m;  /* the values to be read */

    /* Scan through file line by line */
    while (fgets(buf, sizeof(buf), fin) != NULL) {
        lineno++;
        start = lskip(rstrip(buf));     /* chop whites */

        if (*start && *start != '#') {
            /* Remove possible comments at the end */
            end = find_char_or_comment(start, '#');
            if (*end == '#')
                *end = '\0';
            values = rstrip(start);
            /* Not a comment, must be a line containing 6 numbers */
            /* Scan and check if the inputs is correctly formatted */
            if (sscanf(values, "%lf %lf %lf %lf %lf %lf",
                       &t, &d, &q1, &q2, &i, &m) == 6) {
                if (t >= 0 && d > 0 && m >= 0) {
                    *cfginputs =
                        append_pulse_2D(*cfginputs,
                                        new_pulse_2D(N1, N2, t, d,
                                                     -q1, -q2, i, m));
                } else {
                    printf("Durations and time onsets should be positive.\n");
                    printf("Check configuration file for external inputs.\n");
                    error = lineno;
                }
            } else if (!error) {
                /* No '=' found on name=value line  */
                error = lineno;
            }
        }
    }
    fclose(fin);
    return error;
}
