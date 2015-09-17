#ifndef __COLORMAP_H_
#define __COLORMAP_H_ 1

#include <math.h>
/* Data structures for colormap stuff */
struct intval_node {
    double position;
    double value;
};

struct color_desc {
    int n_subintervals;
    struct intval_node *nodes;
};

void interpolate_col(const struct color_desc *color, unsigned int table[][3],
                     const int N_rows, const int col_index);
void fill_cm_lookuptable_jet(unsigned int table[][3], const int N_rows);
void fill_cm_lookuptable_bone(unsigned int table[][3], const int N_rows);

#endif
