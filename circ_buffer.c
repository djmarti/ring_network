#include "circ_buffer.h"

double **allocate_matrix(size_t M, size_t N)
{
    double **a;
    /*  allocate storage for an array of pointers */
    a = malloc(M * sizeof(double *));

    /* for each pointer, allocate storage for an array of doubles */
    for (size_t i = 0; i < M; i++) {
        a[i] = malloc(N * sizeof(double));
    }
    return a;
}

void free_matrix(double **a, size_t M)
{
    for (size_t i = 0; i < M; i++)
        free(a[i]);
    free(a);
}

struct circ_buffer circ_buffer_create(size_t time, size_t Nn)
{
    struct circ_buffer cb;
    cb.size = time + 1;         /* we need one extra slot to diff start from end */
    cb.ncols = Nn;
    cb.image = allocate_matrix(cb.size, Nn);
    cb.start = 0;
    cb.next = 0;
    return cb;
}

void circ_buffer_add_vector(struct circ_buffer *cb, const double *v)
{
    /* Insert row at the approprate place */
    for (size_t j = 0; j < cb->ncols; j++) {
        cb->image[cb->next][j] = v[j];
    }
    cb->next = (cb->next + 1) % cb->size;
    /* Is the buffer full? */
    if (cb->next == cb->start) {
        cb->start = (cb->start + 1) % cb->size;
    }
}

void circ_buffer_free(struct circ_buffer *cb)
{
    free_matrix(cb->image, cb->size);
}
