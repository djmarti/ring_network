#ifndef __CIRC_BUFFER_H_
#define __CIRC_BUFFER_H_ 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* The buffer is basically a matrix with 'size' rows and 'ncols' columns.
 * The vectors of a variable measured at time t are stored as a row in that
 * matrix, and successive vectors are stacked one _under_ the other (as the
 * first-order index increases). That is, time increases as we go to higher row
 * indices.
 */

struct circ_buffer {
    size_t size;                /* Number of rows (time snapshots) */
    size_t ncols;               /* Number of columns (vector component) */
    double **image;             /* The actual chunk of data */
    size_t start;               /* Where the initial index is. Usually this stays at 0 unless
                                   we reach the end of the buffer */
    size_t next;                /* next empty position */
};

double **allocate_matrix(size_t, size_t);
void free_matrix(double **a, size_t M);
struct circ_buffer circ_buffer_create(size_t time, size_t Nn);
void circ_buffer_add_vector(struct circ_buffer *cb, const double *v);
void circ_buffer_free(struct circ_buffer *cb);
#endif
